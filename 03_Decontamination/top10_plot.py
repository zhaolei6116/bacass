import matplotlib.pyplot as plt
import csv
import argparse
import seaborn

def read_tsv_file(filename):
    """
    从TSV文件中读取数据。
    """
    species_names = []
    species_reads = []
    with open(filename, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # 跳过标题行
        for row in reader:
            species_names.append(row[0])
            species_reads.append(int(row[1]))
    return species_names, species_reads

def plot_data(species_names, species_reads):
    """
    使用matplotlib绘制柱状图。
    """
    # plt.style.use('seaborn-darkgrid')
    plt.figure(figsize=(14, 8))
    bars = plt.barh(species_names, species_reads, color='#85C1E9', alpha=0.7)
    plt.ylabel('Species Name')
    plt.xlabel('Species Reads')
    plt.title('Top 10 Species', fontsize=16)
    plt.xticks(fontsize=10)
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 5, bar.get_y() + bar.get_height()/2, f'{width:,}', va='center', ha='left', fontsize=10)
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("top10.png")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot data from a TSV file.')
    parser.add_argument('filename', type=str, help='The path to the TSV file.')
    args = parser.parse_args()
    
    species_names, species_reads = read_tsv_file(args.filename)
    plot_data(species_names, species_reads)

if __name__ == "__main__":
    main()