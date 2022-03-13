"""
=================
Pipeline memory
=================
Increase memory for R scripts

"""

from ruffus import *

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

    
@originate("MM_relapsed_scclassify_model.rds")
def r_scclassify_train(outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript build_scclassify_model.R
    '''

    P.run(statement)

@originate("MM_joint_scclassify_model.rds")
def r_scclassify_train(outfile):
    '''
    Runs Rscript
    '''


    job_memory = "75G"

    statement = '''
    Rscript build_scclassify_model_joint.R
    '''

    P.run(statement)

def main(argv=None):
	if argv is None:
		argv = sys.argv
	P.main(argv)

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))