import pysam
from matplotlib import pyplot as plt
import numpy as np
import argparse

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='',
        description="")

    parser.add_argument("-b", "--bam_file",
                        required=True,
                        help="input bam file")
    parser.add_argument("-p", "--out_prefix",
                        required=True,
                        help="location for output files")
    return parser.parse_args()

args = arg_parser()

bamfile = pysam.AlignmentFile(args.bam_file, "rb")

de_list=[]
for read in bamfile.fetch():

    try:
        tag=str(read.get_tag("de"))
        if float(tag) > 0: 
            de_list.append(tag)
    except:
        print("de tag not found in read")
data=np.array(de_list)
data = data.astype(float)

# write out de list
with open(args.out_prefix + ".de_dist.txt", 'w') as fp:
        for item in de_list:
                    # write each item on a new line
                            fp.write("%s\n" % item)
# Set the bin size
bin_size = 1e-2

# Calculate the number of bins based on the bin size
#num_bins = int(1 / bin_size)
num_bins = int(np.ceil((np.max(data) - np.min(data)) / bin_size))

# Create the histogram
plt.hist(data, bins=num_bins, color = "blue",alpha=1, rwidth=0.5)

plt.xlabel('de Value')
plt.ylabel('Frequency')
plt.title('de dist of reads, binsize 1e-2')

interval = np.arange(0,0.25,0.01) 
plt.set_xticks(xinterval)
plt.set_xticklabels(xinterval)
# Display the plot
plt.savefig(args.out_prefix +".de_dist.png", dpi=600)
