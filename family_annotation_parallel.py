import argparse
import os
import subprocess
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def load_pfam2go(pfam2go_file):
    """Загружает соответствие PFAM-GO из файла pfam2go"""
    pfam2go = defaultdict(list)
    with open(pfam2go_file, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue
            pfam_match = re.search(r'Pfam:(\w+)', line)
            go_terms = re.findall(r'GO:(\d+)', line)
            if pfam_match and go_terms:
                pfam_id = pfam_match.group(1)
                pfam2go[pfam_id].extend([f'GO:{go}' for go in go_terms])
    return pfam2go

def parse_hmmsearch_tblout(tblout_file):
    """Парсит результаты hmmsearch"""
    domains = set()
    with open(tblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 13:
                continue
            pfam_id = parts[3].split('.')[0]
            evalue = float(parts[4])
            if evalue < 1e-5:
                domains.add(pfam_id)
    return domains


def parse_hmmsearch_tblout_besthit(tblout_file):
    """Парсит результаты hmmsearch и возвращает домены с минимальным e-value для каждого белка"""
    best_hits = {}  # {target_id: (pfam_id, model_name, evalue)}
    
    with open(tblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 18:  # Проверяем минимальное количество полей
                continue
            
            # Извлекаем данные из столбцов
            target_id = parts[0]           # Идентификатор белка
            model_name = parts[2]          # Название модели (3-й столбец)
            pfam_id_full = parts[3]        # Полный PFAM ID (например PF00067.24)
            pfam_id = pfam_id_full.split('.')[0]  # Базовый PFAM ID
            evalue = float(parts[4])       # E-value из 7-го столбца
            
            # Обновляем лучший хит для белка
            current = best_hits.get(target_id)
            if not current or evalue < current[2]:
                best_hits[target_id] = (pfam_id, model_name, evalue)
    
    # Считаем частоту доменов
    domain_counts = defaultdict(int)
    for pfam_id, model_name, _ in best_hits.values():
        domain_counts[(pfam_id, model_name)] += 1
    
    # Форматируем результат
    formatted_domains = [
        f"{pfam_id} ({model_name}):{count}" 
        for (pfam_id, model_name), count in domain_counts.items()
    ]
    
    return formatted_domains


def process_family(args, family, pfam2go, output_dir):
    """Обрабатывает одно семейство"""
    family_number = family.replace('family', '')
    fasta_path = os.path.join(args.fasta_dir, f"Family{family_number}_target_genes.fa")
    
    if not os.path.exists(fasta_path):
        return f"{family},File not found,"
    
    # Создаем папку для результатов
    os.makedirs(output_dir, exist_ok=True)
    tblout_file = os.path.join(output_dir, f"{family}_hmmsearch.tbl")
    
    # Запуск hmmsearch с многопоточностью
    cmd = [
        'hmmsearch',
        '--tblout', tblout_file,
        '--noali',
        '--cpu', str(args.threads),
        args.pfam_db,
        fasta_path
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        return f"{family},Error: {e},"
    
    # Парсинг результатов
    domains = parse_hmmsearch_tblout_besthit(tblout_file)
    
    # Сбор GO терминов
    go_terms = set()
    for domain in domains:
        go_terms.update(pfam2go.get(domain.split()[0], []))
    
    return f"{family},{';'.join(domains)},{';'.join(go_terms)}"

def main():
    parser = argparse.ArgumentParser(description='Анализ PFAM и GO с многопоточностью')
    parser.add_argument('--family_list', required=True)
    parser.add_argument('--fasta_dir', required=True)
    parser.add_argument('--pfam_db', required=True)
    parser.add_argument('--pfam2go', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--tmp_dir', default='hmm_results', help='Папка для промежуточных результатов')
    parser.add_argument('--threads', type=int, default=4, help='Количество потоков для hmmsearch')
    args = parser.parse_args()

    # Создаем основную папку для результатов
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Загрузка данных
    with open(args.family_list, 'r') as f:
        families = [line.strip().lower() for line in f]
    
    pfam2go = load_pfam2go(args.pfam2go)

    # Обработка с использованием ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=args.threads) as executor, open(args.output, 'w') as out_f:
        out_f.write("Family,PFAM_Domains,GO_Terms\n")
        futures = []
        for family in families:
            futures.append(
                executor.submit(
                    process_family,
                    args,
                    family,
                    pfam2go,
                    args.tmp_dir
                )
            )
        
        for future in futures:
            result = future.result()
            out_f.write(result + "\n")

if __name__ == "__main__":
    main()
