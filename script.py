import urllib.request
import gzip
import shutil
import os
import subprocess
import csv

# Download and unzip the file
url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"

downloaded_file = "clinvar.vcf.gz"
unzipped_file = "./snpEff/data/clinvar.vcf"

if not os.path.exists(unzipped_file):
    # Download the file
    print("Downloading file from:", url)
    urllib.request.urlretrieve(url, downloaded_file)
    print("Download complete.")

    print("Unzipping", downloaded_file)
    with gzip.open(downloaded_file, 'rb') as f_in:
        with open(unzipped_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("Unzipping complete.")
    os.remove(downloaded_file)
    print("Data downloaded and unzipped successfully.")
else:
    print("File already exists. Skipping download.")

#Annotate with snpEff
bash_command = "java -Xmx8g -jar ./snpEff/snpEff.jar -v -c ./snpEff/snpEff.config GRCh37.75 ./snpEff/data/clinvar.vcf > ./snpEff/data/clinvar_annotated.vcf"
print("Running command:", bash_command)
os.system(bash_command)
print("Annotation complete.")


data = []
with open("./snpEff/data/clinvar_annotated.vcf", 'r') as infile:
    for line in infile:
        if not line.startswith('#'):  # Skip header lines
            data.append(line.strip().split('\t'))
            
annotation = [[],[],[]]      
chromosome = [[],[]]
disease = []
variant = [[],[],[]]      
interpretation = [[],[],[],[],[]]
hgvs = []
locationInfo = [[],[]]
            
for row in data:
    chromosome[0].append(row[0])
    chromosome[1].append("GRCh37.75")
    variant[1].append(row[4])
    locationInfo[0].append(row[1])
    locationInfo[1].append(row[3])
    for i in row[-1].split(';'):
        if i.startswith('CLNDN'):
            disease.append(i.split('=')[1])
        if i.startswith('RS'):
            variant[0].append(i.split('=')[1])
        if i.startswith('CLNVC'):
            variant[2].append(i.split('=')[1])
        if i.startswith('CLINSIG'):
            interpretation[0].append(i.split('=')[1])
        if i.startswith('CLNVI'):
            interpretation[1].append(i.split('=')[1])
        if i.startswith('ORIGIN'):
            interpretation[2].append(i.split('=')[1])
        if i.startswith('CLNREVSTAT'):
            interpretation[3].append(i.split('=')[1])
        if i.startswith('CLNDISDB'):
            interpretation[4].append(i.split('=')[1])
        if i.startswith('CLNHGVS'):
            hgvs.append(i.split('=')[1])
        if i.startswith('ANN'):
            annotation[0].append(i.split('|')[2])
            annotation[1].append(i.split('|')[-1])
            annotation[2].append(i.split('|')[0][4:])
        

with open('./snpEff/data/chromosome.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['chromosome', 'assembly'])
    for i in range(len(chromosome[0])):
        clinvar_writer.writerow([chromosome[0][i], chromosome[1][i]])

with open('./snpEff/data/variant.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['variant_rs_id', 'alt', 'variant_type'])
    for i in range(len(variant[0])):
        clinvar_writer.writerow([variant[0][i], variant[1][i], variant[2][i]])

with open('./snpEff/data/locationInfo.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['pos', 'ref'])
    for i in range(len(locationInfo[0])):
        clinvar_writer.writerow([locationInfo[0][i], locationInfo[1][i]])

with open('./snpEff/data/disease.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['preferred_name'])
    for i in range(len(disease)):
        clinvar_writer.writerow([disease[i]])

with open('./snpEff/data/interpretation.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['clinical_significance', 'method', 'variant_origin', 'review_status', 'submitter'])
    for i in range(len(interpretation[0])):
        clinvar_writer.writerow([interpretation[0][i], interpretation[1][i], interpretation[2][i], interpretation[3][i], interpretation[4][i]])
        
with open('./snpEff/data/hgvs.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['hgvs'])
    for i in range(len(hgvs)):
        clinvar_writer.writerow([hgvs[i]])

with open('./snpEff/data/annotation.csv', mode='w') as clinvar_file:
    clinvar_writer = csv.writer(clinvar_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    clinvar_writer.writerow(['impact', 'consequence', 'allele'])
    for i in range(len(annotation[0])):
        clinvar_writer.writerow([annotation[0][i], annotation[1][i], annotation[2][i]])