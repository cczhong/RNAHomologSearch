import random
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write


def get_consensus_from_stockholm(msa_file):
    consensus = None
    with open(msa_file, "r") as f:
        for line in f:
            if line.startswith("#=GC RF"):
                consensus = line.strip().split(maxsplit=2)[2]
                break
    if not consensus:
        raise ValueError("Consensus sequence not found in the MSA file.")
    return consensus


def generate_challenge_sequences(msa_file, file_format, num_sequences, seq_length=100, min_align_length=40):
    try:
        alignment = AlignIO.read(msa_file, file_format)
        consensus = get_consensus_from_stockholm(msa_file)

        base_name = os.path.splitext(os.path.basename(msa_file))[0]
        output_fasta = f"{base_name}_output.fasta"
        output_csv = f"{base_name}_output.csv"

        fasta_records = []
        csv_records = []

        for i in range(num_sequences):
            selected_record = random.choice(alignment)
            sequence_id = selected_record.id
            sequence = str(selected_record.seq)

            valid_positions = [j for j, c in enumerate(consensus) if c != "-" and sequence[j] != "-"]
            if len(valid_positions) < min_align_length:
                raise ValueError("No valid subsequences meeting the minimum align length.")

            start_index = random.choice(valid_positions)
            end_index = start_index + min_align_length

            fragment = sequence[start_index:end_index]
            if len(fragment) < min_align_length:
                continue

            total_pad = seq_length - len(fragment)
            left_pad = "".join(random.choices("AUCG", k=total_pad // 2))
            right_pad = "".join(random.choices("AUCG", k=total_pad - len(left_pad)))

            challenge_seq = left_pad + fragment + right_pad

            fasta_records.append(SeqRecord(Seq(challenge_seq), id=f"Challenge_{i}", description=""))
            csv_records.append((sequence_id, start_index, end_index, challenge_seq))

        with open(output_fasta, "w") as fasta_out:
            write(fasta_records, fasta_out, "fasta")

        with open(output_csv, "w") as csv_out:
            csv_out.write("Sequence_ID,Start,End,Challenge_Sequence\n")
            for rec in csv_records:
                csv_out.write(",".join(map(str, rec)) + "\n")

        print(f"Challenge sequences written to {output_fasta} and boundaries to {output_csv}.")

    except Exception as e:
        print(f"Error: {e}")

def main():
    file_name = input("Enter the file name: ")
    generate_challenge_sequences(file_name, "stockholm", 5, 100, 40)

main()