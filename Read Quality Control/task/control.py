# write your code here
from collections import Counter
from statistics import mean
import gzip


# from pathlib import Path

# root = Path(__file__).parent.parent.parent
# path_ = f"{root}/{user_input}"

def calculate_quality(path_):
    seq_set = set()
    repeats = list()
    reads_with_ns = list()
    all_ns = list()

    with gzip.open(path_, "rb") as file:
        data = file.read().decode("utf-8").split("\n")

    count_of_length = dict()
    gc_content = list()

    def calculate_content_value(seq):
        counts = dict(Counter(seq))
        gc_value = sum(v for k, v in counts.items() if k in ["G", "C"])
        all_ = sum(v for k, v in counts.items())
        Ns = sum(v for k, v in counts.items() if k == "N") / all_ * 100
        reads_with_ns.append(1) if Ns > 0.5 else None
        all_ns.append(Ns)

        return gc_value / all_ * 100

    for i in range(1, len(data), 4):

        sequence = data[i]

        seq_set.add(sequence)

        gc_content.append(calculate_content_value(sequence))

        length = len(sequence)

        if length not in count_of_length:
            count_of_length.update({length: 1})
        else:
            count_of_length.update({length: count_of_length[length] + 1})

    sorted_count_of_length = sorted(count_of_length.items(), key=lambda x: x[0])

    n_reads = sum([v for k, v in sorted_count_of_length])
    average_reads = sum([k * v for k, v in sorted_count_of_length]) / n_reads

    return {
        "n_reads": n_reads,
        "average_reads": average_reads,
        "seq_set": seq_set,
        "reads_with_ns": reads_with_ns,
        "gc_content": gc_content,
        "all_ns": all_ns,
    }


paths = list()

for _ in range(3):
    paths.append(input())

responses = list()

for path in paths:
    responses.append(calculate_quality(path))

best_response = sorted(responses, key=lambda x: x["gc_content"])[-1]

print(f"Reads in the file = {best_response['n_reads']}:")
print(f"Reads sequence average length = {round(best_response['average_reads'], 0)}")
print(f"Repeats = {best_response['n_reads'] - len(best_response['seq_set'])}")
print(f"Reads with Ns = {sum(best_response['reads_with_ns'])}")
print(f"GC content average = {round(mean(best_response['gc_content']), 2)}%")
print(f"Ns per read sequence = {round(sum(best_response['all_ns']) / best_response['n_reads'], 2)}%")
