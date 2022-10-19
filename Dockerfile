FROM ashedpotatoes/sranwrp:1.0.3
ADD cut-megafile-sept.py ./scripts/
ENV PATH=/scripts:/root/miniconda3/bin:/root/miniconda3/condabin:/root/miniconda3/bin:/bin:/root/edirect:/sra-tools-3.0.0:/ncbi-vdb-3.0.0:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin