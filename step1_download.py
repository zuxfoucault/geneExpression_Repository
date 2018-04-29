import os
import time
from subprocess import call

# allen brain institute dataset: 6 donors names and links
donors_name=['normalized_microarray_donor9861',
             'normalized_microarray_donor10021',
             'normalized_microarray_donor12876',
             'normalized_microarray_donor14380',
             'normalized_microarray_donor15496',
             'normalized_microarray_donor15697']

donors_links=['http://human.brain-map.org/api/v2/well_known_file_download/178238387',
              'http://human.brain-map.org/api/v2/well_known_file_download/178238373',
              'http://human.brain-map.org/api/v2/well_known_file_download/178238359',
              'http://human.brain-map.org/api/v2/well_known_file_download/178238316',
              'http://human.brain-map.org/api/v2/well_known_file_download/178238266',
              'http://human.brain-map.org/api/v2/well_known_file_download/178236545']


donors_links_T1=['http://human.brain-map.org/api/v2/well_known_file_download/157722636',
                'http://human.brain-map.org/api/v2/well_known_file_download/157723301',
                'http://human.brain-map.org/api/v2/well_known_file_download/157722290',
                'http://human.brain-map.org/api/v2/well_known_file_download/157721937',
                'http://human.brain-map.org/api/v2/well_known_file_download/162021642',
                'http://human.brain-map.org/api/v2/well_known_file_download/157682966']

donors_links_FS=['https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor9861.zip?sequence=1&isallowed=y',
                 'https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor10021.zip?sequence=2&isallowed=y',
                 'https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor12876.zip?sequence=3&isallowed=y',
                 'https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor14380.zip?sequence=4&isallowed=y',
                 'https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor15496.zip?sequence=5&isallowed=y',
                 'https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/donor15697.zip?sequence=6&isallowed=y']

download_folder = 'AIBS_map/Allen_FS/'
if not os.path.exists(download_folder):
    os.makedirs(download_folder)
# download table with structures names
link_table="https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/structurelist_rrgedit.csv?sequence=7&isallowed=y"
#call(["wget", "-O", "AIBS_map/Allen_FS/StructureList_RRGedit.csv", link_table])


# download fsaverage subject
link_fsaverage="https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/fsaveragesubp.zip?sequence=9&isallowed=y"
#call(["wget", "-O", "AIBS_map/Allen_FS/fsaverageSubP.zip", link_fsaverage])
#time.sleep(15)
#call(["unzip", "-o", "AIBS_map/Allen_FS/fsaverageSubP.zip", "-d", "AIBS_map/Allen_FS/fsaverageSubP"])

# Add auxiliar folder to include required functions
# specific for matlab
#matlabpath_auxiliar=which('main_sampleMatching');
#matlabpath_auxiliar=[path_auxiliar(1:end-22) '/auxiliar/'];
#matlabaddpath(path_auxiliar);

# Loop across donors
download_folder = "AIBS_map/downloaded/"
if not os.path.exists(download_folder):
    os.makedirs(download_folder)
for ifol in range(len(donors_name)):
    donor_name=donors_name[ifol]
    donor_link=donors_links[ifol]
    donor_link_T1=donors_links_T1[ifol]
    donor_link_FS=donors_links_FS[ifol]

    # Download freesurfer folder of each donor
    call(["wget", "-O", "AIBS_map/Allen_FS/"+donor_name+".zip", donor_link_FS])
    time.sleep(15)
    call(["unzip", "-o", "AIBS_map/Allen_FS/"+donor_name+".zip", "-d", "AIBS_map/Allen_FS/"])

    # Download data from AIBS web
    call(["wget", "-O", "AIBS_map/downloaded/"+donor_name+".zip", donor_link])
    time.sleep(15)
    call(["unzip", "-o", "AIBS_map/downloaded/"+donor_name+".zip", "-d", "AIBS_map/downloaded/"+donor_name])

    call(["wget", "-O", "AIBS_map/downloaded/"+donor_name+"/T1.nii.gz", donor_link_T1])
    time.sleep(15)

    filename="AIBS_map/downloaded/" + donor_name + "/Probes.csv"

#    # Import Probe file
#    [probe_id,probe_name,gene_id,gene_symbol,gene_name,entrez_id,chromosome]=import_probe(filename);
#    gene_symbol_orig=gene_symbol;
#    %Reannot probes using Richiardi et al 2015 table
#    reannot=importdata('auxiliar/Richiardi_Data_File_S2_Extended.csv');
#    reannot_probe=reannot.textdata(2:end,1);
#    reannot_id=reannot.data;
#    reannot_genes=reannot.textdata(2:end,4);
#    for in=1:numel(probe_name)
#        reannot_matching=find(strcmp(probe_name{in},reannot_probe));
#        if (not(isempty(reannot_matching))) && not(isempty(reannot_genes{reannot_matching}))
#            gene_symbol{in}=reannot_genes{reannot_matching};
#            %entrez_id(in)=reannot_id(reannot_matching);
#        end
#    end
#    %Do not include Probes without entrez_id
#    find_complete_valid_cond1=find(isnan(entrez_id)==0);
#
#    %Do not include Probes not associated to genes id
#    find_complete_valid_cond2=find(not(strcmp(gene_symbol,'na')));
#    find_complete_valid=intersect(find_complete_valid_cond1,find_complete_valid_cond2);
#    gene_symbol_valid=gene_symbol(find_complete_valid);
#    gene_valid_notRep=unique(gene_symbol_valid);
#
#    %Total number of genes should be 20737 (This was true before reannotation)
#    %Total number of genes should be 20647 (This was true before reannotation)
#
#    if numel(gene_valid_notRep)~=20647
#        error('Wrong number of genes');
#    end
#
#    %Extract valid probes
#    res_probid_val=probe_id(find_complete_valid);
#    res_gene_symbol_val=gene_symbol(find_complete_valid);
#    res_u_gene_symbol=unique(res_gene_symbol_val);
#
#    %Look for all the probes matching gene label (include only valid genes)
#    for il=1:numel(res_gene_symbol_val)
#        pos=find(strcmp(res_gene_symbol_val{il},res_u_gene_symbol));
#        res_gene_symbol_val_tonum(il)=pos;
#    end
#
#    %Look for all the probes matching gene label
#    for il=1:numel(gene_symbol)
#        pos=find(strcmp(gene_symbol{il},res_u_gene_symbol));
#        if isempty(pos)
#            pos=-1;
#        end
#        res_gene_symbol_tonum(il)=pos;
#    end
#
#    %Save output
#    path_output=['AIBS_map/downloaded/' donor_name '/probe2gene/'];
#    [s1 s2 s3]=mkdir(path_output);
#
#    %List that tranform from Probe to gene (only valid)
#    save([path_output 'res_gene_symbol_val_tonum.mat'],'res_gene_symbol_val_tonum');
#
#    %List of all valid probes
#    save([path_output 'res_probid_val.mat'],'res_probid_val');
#
#    %List that tranform from Probe to gene
#    save([path_output 'res_gene_symbol_tonum.mat'],'res_gene_symbol_tonum');
#
#    %List of all probes
#    save([path_output 'probe_id.mat'],'probe_id');
#
#    %List of all valid genes name (unique)
#    save([path_output 'res_u_gene_symbol.mat'],'res_u_gene_symbol');
#
#    display(['Subject ' donor_name ' Done!']);
#
#end
#display(['Probe to gene completed']);
#end
