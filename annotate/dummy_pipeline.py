"""
=================
Pipeline memory
=================
Increase memory for R scripts

"""

from ruffus import *

import os
import sys
import cgatcore.pipeline as P
import cgatcore.experiment as E



# Load options from the config file

@originate("GSM5687372_filtered_SeuratObject.rds")
def r_gsm56(outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript create_seurat_GSM56.R
    '''

    P.run(statement)

filtered_suffixes = "*_integrated_SeuratObject.rds"
SEURAT_OBJECTS = tuple([os.path.join("RDS_objects.dir",filtered_suffixes)])

@transform(SEURAT_OBJECTS, 
          regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
          r"RDS_objects.dir/\1_annotated_MM_naive_classifier_SeuratObject.rds")
def r_scclassify_MM_naive(infile,outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript annotate_MM_naive_classifier.R -i %(infile)s -o %(outfile)s
    '''

    P.run(statement)

@transform(SEURAT_OBJECTS, 
          regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
          r"RDS_objects.dir/\1_annotated_MM_relapse_classifier_SeuratObject.rds")
def r_scclassify_MM_relapse(infile,outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript annotate_MM_relapse_classifier.R -i %(infile)s -o %(outfile)s
    '''

    P.run(statement)

@transform(SEURAT_OBJECTS, 
        regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
        r"RDS_objects.dir/\1_annotated_MM_joint_classifier_SeuratObject.rds")
def r_scclassify_MM_joint(infile,outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript annotate_MM_joint_classifier.R -i %(infile)s -o %(outfile)s
    '''

    P.run(statement)
    
    
@transform(SEURAT_OBJECTS, 
        regex("RDS_objects.dir/(\S+)_integrated_SeuratObject.rds"),
        r"RDS_objects.dir/\1_clustifyr_annotated_MM_classifier_SeuratObject.rds")
def r_clustifyr_MM_joint(infile,outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript clustifyr_MM_annotate.R -i %(infile)s 
    '''

    P.run(statement)
    
def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))