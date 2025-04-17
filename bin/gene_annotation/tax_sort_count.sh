#!/bin/bash

# 检查是否提供了必要的参数
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <sorted_output_file> <final_output_file>"
    exit 1
fi

# 输入文件名
input_file="$1"

# 排序后的中间文件名（由外接设定）
sorted_output_file="$2"

# 最终输出文件名（包含前5物种和 "Other" 的统计）
final_output_file="$3"

# 使用 awk 统计每个物种的 count 总数，并按顺序排序
awk -F'\t' '
 NR > 1  {  # 忽略表头
    species_count[$NF] += $2  # 按最后一列（species）累加 count
}
END {
    for (species in species_count) {
        print species "\t" species_count[species]  # 使用 Tab 分隔输出
    }
}' "$input_file" | sort -t$'\t' -k2,2nr -k1,1 > "$sorted_output_file"

# 处理排序后的结果，提取前5个物种并将剩余归类为 "Other"
{
    echo -e "Species\tCount" > "$final_output_file"  # 写入头部信息
    i=0
    total_other=0
    while IFS=$'\t' read -r species count; do
        if [ $i -lt 5 ]; then
            echo -e "$species\t$count" >> "$final_output_file"
        else
            total_other=$((total_other + count))
        fi
        ((i++))
    done
} < "$sorted_output_file"

# 将 "Other" 类别的统计数据追加到输出文件
echo -e "Other\t$total_other" >> "$final_output_file"

# 清理临时文件（可选）
# rm "$sorted_output_file"

echo "Processing complete. Results saved to '$final_output_file'."
