import pandas as pd
import fileinput


def convert_to_family_file(ortho_result):
    """Функция переводит реультаты ProyeinOrtho
    в входную таблицу CAFE"""
    ortho_result = pd.read_csv(ortho_result, sep='\t').iloc[:, 3:]
    weight_result = ortho_result.map(lambda x: 0 if x == '*' else len(x.split(',')))

    list_family_ids = ["family" + str(i) for i in range(1, len(ortho_result)+1)]

    weight_result.insert(0, "Desc", '(null)')
    weight_result.insert(1, "Family ID", list_family_ids )
    return weight_result


def mafft_prep(ortho_result):
    """Вычленяет из результатов ProteinOrtho только строки,
    где 1 ортолог у каждого вида"""
    ortho_result = pd.read_csv(ortho_result, sep='\t')
    num = ortho_result.iloc[:, 3:].shape[1]
    ortho_result = ortho_result.loc[lambda df: df['Genes'] == num]
    ortho_result = ortho_result[ortho_result.ne('*').all(axis=1)]

    return ortho_result.sort_values(by='Alg.-Conn.', ascending=False).drop(columns=['Genes', 'Alg.-Conn.', '# Species'])

def input_mafft(one_copy_table, row):
    """Функция печатает id генов чтобы отправить их в samtools,
    возвращает список кортежей с парами ген-вид для замены"""
    prep_table = one_copy_table.iloc[[row]]
    for i in range(len(prep_table)):
        print(prep_table.iloc[i])
    with open(f'{prep_table.index[0]}fam_proteins.txt', 'w', encoding='utf-8') as file:
        for item in prep_table.iloc[0].tolist():
            file.write(f"{item}\t")

    return list(one_copy_table.iloc[row, :].items())

def replace_names(mafft_res, replacement_list):
    """Функция заменяет id гена на название вида"""
    for item in replacement_list:
        new_line = '>' + str(item[0])
        old_line = '>' + str(item[1])
        for line in fileinput.input(mafft_res, inplace=True):
            print(line.replace(old_line, new_line), end='')

def significant(cafe_sign, ortho_res):
    """Функция получает на вход файл со значительно изменившимися семействами из
    CAFE и таблицу результатов Proteinortho. Возвращает папку с текстовыми
    файлами, где в каждом файле лежит набор названий генов, разделенных табуляцией """
    cafe_sign = pd.read_csv(cafe_sign, sep = '\t', header = None)
    fam_list = cafe_sign.iloc[:, 0].tolist()

    for x in range(len(fam_list)):
        fam_list[x] = int(fam_list[x].split('y')[1]) - 1
    ortho_res = pd.read_csv(ortho_res, sep = '\t')
    sign_ortho = ortho_res.loc[fam_list].iloc[:, 3:]
    print(sign_ortho)
    for i, ind in zip(range(len(sign_ortho)), fam_list):
        row_df = sign_ortho.iloc[[i]]
        genes_list = []
        for j in range(row_df.shape[1]):
            if row_df.iloc[0, j] != '*':
                genes_list += row_df.iloc[0, j].split(',')
        with open(f'./significant_fam_final/Family{ind + 1}_genes.txt', 'w') as f:
            print(row_df.index)
            for value in genes_list:
                f.write(f"{value}\t")

def incr_decr(base_count, base_change, cafe_sign):
    """Функция получает на вход файлы CAFE, на выходе определяет
    судьбу семейства генов для каждого таксона"""
    # Загрузка файлов
    count_df = pd.read_csv(base_count, sep='\t')
    change_df = pd.read_csv(base_change, sep='\t')

    # Настроим FamilyID как индекс для удобства
    count_df.set_index('FamilyID', inplace=True)
    change_df.set_index('FamilyID', inplace=True)

    cafe_sign = pd.read_csv(cafe_sign, sep = '\t', header = None)
    count_df = count_df.loc[cafe_sign.iloc[:, 0]]
    change_df = change_df.loc[cafe_sign.iloc[:, 0]]

    # Списки для хранения результатов
    disappeared = []
    appeared = []
    increased = []
    decreased = []
    micro_25 = set()

    # Перебираем семейства
    for family in count_df.index:
        counts = count_df.loc[family]
        changes = change_df.loc[family]

        for taxon, (count, change) in zip(counts.index, zip(counts.values, changes.values)):
            if count == 0 and change < 0:
                disappeared.append((family, taxon, change))
            if count > 0 and change > 0 and count == change:
                appeared.append((family, taxon, change))
            if count > 0 and change > 0:
                increased.append((family, taxon, change))
            if count > 0 and change < 0:
                decreased.append((family, taxon, change))

    # Вывод результатов
    print("\nСемейства, которые исчезли:")
    for fam, tax, chg in disappeared:
        if tax == '<25>':
            print(f"{fam}, Таксон: Microsporidia, Исчезло генов: {chg}")
            micro_25.add(fam)

    #print("\nСемейства, которые появились:")
    for fam, tax, chg in appeared:
        print(f"{fam}, Таксон: {tax}, приобретено генов: {chg}")

    print("\nСемейства, которые рисширились:")
    for fam, tax, chg in increased:
        if tax == '<25>':
            print(f"{fam}, Таксон: {tax}, Приобретено генов : {chg}")
            micro_25.add(fam)

    print("\nСемейства, которые уменьшились:")
    for fam, tax, chg in decreased:
        if tax == '<25>':
            print(f"{fam}, Tаксон: Microsporidia, Потеряно генов: {chg}")
            micro_25.add(fam)

    def extract_number(s):
        return int(''.join(filter(str.isdigit, s)))

    with open('25_taxa.txt', 'w', encoding='utf-8') as file:
        for item in sorted(micro_25, key=extract_number):
            file.write(f"{item}\n")

