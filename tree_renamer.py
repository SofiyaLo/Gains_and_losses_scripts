import re

def load_name_mapping(mapping_file):
    mapping = {}
    with open(mapping_file, 'r', encoding='utf-8') as f:
        for line in f:
            old, new = line.strip().split('\t')
            mapping[old] = new
    return mapping

def escape_newick_name(name):
    # Экранируем специальные символы для точного поиска
    return re.escape(name)

def rename_tree_labels_precise(nwk_file, mapping, output_file):
    with open(nwk_file, 'r', encoding='utf-8') as f:
        tree = f.read()

    # Сортируем по убыванию длины, чтобы избежать частичных замен
    sorted_names = sorted(mapping.keys(), key=len, reverse=True)

    for old_name in sorted_names:
        escaped_old = escape_newick_name(old_name)
        new_name = mapping[old_name]

        # Заменяем только полные имена (перед `:`, `,`, `)`)
        pattern = rf'(?<=\(|,){escaped_old}(?=[:),])'
        tree = re.sub(pattern, new_name, tree)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(tree)

    print("Замена завершена. Сохранено в", output_file)

# Пример использования
if __name__ == "__main__":
    rename_tree_labels_precise(
        nwk_file=r'C:\Users\Professional\Downloads\dip_final.nwk',
        mapping=load_name_mapping(r'C:\Users\Professional\Desktop\Таблица_видов_на_русский.txt'),
        output_file="tree_renamed.nwk"
    )