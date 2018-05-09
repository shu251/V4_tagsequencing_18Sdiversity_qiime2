import sys
from subprocess import call

manifest = sys.argv[1]
try:
  primer = sys.argv[2]
except:
  primer = '/usr/local/bioinf/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa'
minlen = '100'

forward = dict()
reverse = dict()
with open(manifest) as fh:
  for line in fh:
    if line.strip().split(',')[2] == 'forward':
      forward[line.strip().split(',')[0]] = line.strip().split(',')[1]
    else:
      reverse[line.strip().split(',')[0]] = line.strip().split(',')[1]

call(["mkdir", "trimmed_seq_dir"])
output = open("trimmed_manifest.txt", "w")
log = open("trimming.log", "w")

output.write("sample-id,absolute-filepath,direction\n")
for sample in forward:
  output.write(sample+','+'$PWD/trimmed_seq_dir/'+sample+'_R1_PE.fastq.gz'+','+'forward'+"\n")
  output.write(sample+','+'$PWD/trimmed_seq_dir/'+sample+'_R2_PE.fastq.gz'+','+'reverse'+"\n")
  call(["java", "-jar", "/usr/local/bioinf/Trimmomatic-0.32/trimmomatic-0.32.jar", "PE", forward[sample], reverse[sample], 'trimmed_seq_dir/'+sample+'_R1_PE.fastq.gz', 'trimmed_seq_dir/'+sample+'_R1_SE.fastq.gz', 'trimmed_seq_dir/'+sample+'_R2_PE.fastq.gz', 'trimmed_seq_dir/'+sample+'_R2_SE.fastq.gz', "LEADING:10", "TRAILING:10", "SLIDINGWINDOW:10:25", "MINLEN:"+minlen, "ILLUMINACLIP:"+primer+":2:30:10"], stderr = log)

output.close()
log.close()
