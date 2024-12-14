#Install cell ranger
cd final_project
export PATH="/Program/cellranger-6.0.0:$PATH"
cellranger
Perform a sitecheck
cellranger sitecheck > sitecheck.txt
cellranger testrun --id=check_install

#Using cell ranger for assembly and gene counting 
cellranger count --id=p40xqy \
--fastqs=/Program/data/GS \
--transcriptome= /Program/reference/reference_mouse/refdata-gex-mm10-2020-A \
--jobmode=local \
--include-introns \
--localcores= 36

ls -1 p40xqy/outs

