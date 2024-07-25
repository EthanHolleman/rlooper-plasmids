from Bio import SeqIO
from pathlib import Path

all_vr_seqs_path = '../resources/vr_seqs/complete_inserts.fa'
OUTPUT_DIR = '../resources/sequences'

def extract_vr_name(record):
    return record.description.split(' ')[-1].split('_')[-1]


def make_rlooper_safe_header(record):
    
    length = len(record.seq)
    vr_name = extract_vr_name(record)
    header = f"{vr_name} range=VRFIXED:1-{length} 5'pad=0 3'pad=0 strand=+ repeatMasking=none"
    
    return header


def write_updated_record(record, vr_name):

    filepath = Path(OUTPUT_DIR).joinpath(vr_name).with_suffix('.fa')
    print(filepath)
    SeqIO.write([record], str(filepath), 'fasta')

def process_record(record):
    
    vr_name = extract_vr_name(record)
    record.id = ''
    record.description = make_rlooper_safe_header(record)
    write_updated_record(record, vr_name)
    


def main():
    
    records = SeqIO.parse(all_vr_seqs_path, 'fasta')
    for each_record in records:
        process_record(each_record)



if __name__ == '__main__':
    main()