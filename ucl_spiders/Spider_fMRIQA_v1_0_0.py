'''
    Author:         Benjamin Yvernault
    contact:        b.yvernault@ucl.ac.uk
    Spider name:    fMRIQA
    Spider version: 1.0.0
    Creation date:  2015-10-07 11:17:50.588519
    Purpose:        Execute a script on fMRI data to generate quality assurance
'''

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Execute a script on fMRI data to generate quality assurance"
__spider_name__ = "fMRIQA"
__version__ = "1.0.0"
__modifications__ = "2015-10-07 11:17:50.588519 - Original write"

# Python packages import
import os
import sys
from dax import XnatUtils, spiders, RESULTS_DIR, ScanSpider

def parse_args():
    '''
    Argument parser for the spider input variables
    by default (set in get_session_argparser):
        -p       : proj_label
        -s       : subj_label
        -e       : sess_label
        -c       : scan_label
        -d       : temp_dir
        --suffix : suffix (suffix for assessor proctype on XNAT)
        --host : host for XNAT (default: XNAT_HOST env variable)
        --user : user on XNAT (default: XNAT_USER env variable)
    your arguments:
        --spm
        --fmatlab

    :return: argument parser object created by parse_args()
    '''
    ap = spiders.get_scan_argparser("fMRIQA", "Execute a script on fMRI data to generate quality assurance")
    ap.add_argument("--spm", dest="spm", help="Spm path.", required=True)
    ap.add_argument("--fmatlab", dest="fmatlab", help="Matlab folder containing the needed code.", required=True)
    return ap.parse_args()

class Spider_fMRIQA(ScanSpider):
    '''
        Scan Spider class to do: Execute a script on fMRI data to generate quality assurance

        :param spider_path: spider file path
        :param jobdir: directory for temporary files
        :param xnat_project: project ID on XNAT
        :param xnat_subject: subject label on XNAT
        :param xnat_session: experiment label on XNAT
        :param xnat_scan: scan label on XNAT
        :param xnat_host: host for XNAT if not set in environment variables
        :param xnat_user: user for XNAT if not set in environment variables
        :param xnat_pass: password for XNAT if not set in environment variables
        :param suffix: suffix to the assessor creation
    '''
    def __init__(self, spider_path, jobdir, spm, fmatlab, xnat_project, xnat_subject, xnat_session, xnat_scan,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        super(Spider_fMRIQA, self).__init__(spider_path, jobdir, xnat_project, xnat_subject, xnat_session, xnat_scan,
                                            xnat_host, xnat_user, xnat_pass, suffix)
        self.inputs = list()
        self.spm = spm
        self.fmatlab = fmatlab

    def pre_run(self):
        '''
            Method to download data from XNAT

            :param argument_parse: argument parser object return by parse_args()
        '''
        resource = 'NIFTI' #resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        self.inputs = self.download(self.xnat_scan, resource, folder)
        if not self.inputs:
            raise "No NIFTI files downloaded from XNAT for the scan."
        else:
            for index, image_path in enumerate(self.inputs):
                if image_path.endswith(".nii.gz"):
                    os.system("gzip -f -d %s" % (image_path))
                    self.inputs[index] = image_path[:-3]

    def run(self):
        '''
            Method running the process for the spider on the inputs data
        '''
        matlabdir = os.path.join(self.jobdir, 'MATLAB/')
        if not os.path.exists(matlabdir):
            os.makedirs(matlabdir)
        outputdir = os.path.join(self.jobdir, 'Outputs/')
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)

        matlab_script = os.path.join(matlabdir, 'callfMRIQA_v2.m')
        f = open(matlab_script, "w")
        try:
            lines=[ '% Matlab Script to call vufMRIQAGUI function\n',
                    'addpath(genpath(\''+str(self.fmatlab)+'\'));\n',
                    'outputdir = \''+str(outputdir)+'\';\n',
                    'imgfile=\''+str(self.inputs[0])+'\';\n',
                    'spm_path=\''+str(self.spm)+'\';\n',
                    'project=\''+str(self.xnat_project)+'\';\n',
                    'subject=\''+str(self.xnat_subject)+'\';\n',
                    'session=\''+str(self.xnat_session)+'\';\n',
                    'scan=\''+str(self.xnat_scan)+'\';\n',
                    'fMRIQA_Pipeline_v2(imgfile,outputdir,spm_path,project,subject,session,scan);\n'
                  ]
            f.writelines(lines)
        finally:
            f.close()

        #Running Matlab script:
        XnatUtils.run_matlab(matlab_script, True)
        print '===================================================================\n'

    def finish(self):
        '''
            Method to copy the results in the Spider Results folder dax.RESULTS_DIR
        '''
        self.time_writer('Results saved in folder: %s' % (RESULTS_DIR))
        pdf_name = "%s-x-%s-x-%s-x-%s-x-fMRIQA_v2.ps" % (self.xnat_project, self.xnat_subject, self.xnat_session, self.xnat_scan)
        pdf_fname = os.path.join(self.jobdir,'Outputs',pdf_name)
        os.rename(os.path.join(self.jobdir,'Outputs','fMRIQA_v2.ps'),pdf_fname)
        results_dict = {'PDF': pdf_fname,
                        'STATS': os.path.join(self.jobdir,'Outputs','fMRIQA_v2_stats.txt'),
                        'MATLAB': os.path.join(self.jobdir,'MATLAB')
                        }
        self.upload_dict(results_dict)
        self.end()

if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_fMRIQA(spider_path=sys.argv[0],
                               jobdir=args.temp_dir,
                               spm=args.spm,
                               fmatlab=args.fmatlab,
                               xnat_project=args.proj_label,
                               xnat_subject=args.subj_label,
                               xnat_session=args.sess_label,
                               xnat_scan=args.scan_label,
                               xnat_host=args.host,
                               xnat_user=args.user,
                               xnat_pass=None,
                               suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    spider_obj.finish()
