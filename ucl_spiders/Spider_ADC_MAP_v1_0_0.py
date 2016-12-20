"""Spider_ADC_MAP.

Author:         Benjamin Yvernault & Vasia Papoutsaki
contact:        b.yvernault@ucl.ac.uk
Spider name:    ADC_MAP
Spider version: 1.0.0
Creation date:  2016-12-01 13:49:23.748867
Purpose:        Compute ADC Maps for HandN project
"""

# Python packages import
import os
import sys
import shutil
import subprocess as sb
from dax import XnatUtils, spiders, ScanSpider

__author__ = "Benjamin Yvernault & Vasia Papoutsaki"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Compute ADC Maps for HandN project"
__spider_name__ = "ADC_MAP"
__version__ = "1.0.0"
__modifications__ = """2016-12-01 13:49:23.748867 - Original write"""

DEFAULT_ADC_TEMPLATE = """
addpath(genpath('{matlab_code}'));
XNAT_ADCsMaps('{input_path}',\
'{output}',\
'{pdf_name}');
"""

FSLSWAP_VAL = {0: 'x',
               1: 'y',
               2: 'z'}


def parse_args():
    """Argument parser for the spider input variables.

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
        ...

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_scan_argparser("ADC_MAP", __purpose__)
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help="Matlab code folder where is \
launch_AMICO_for_INNOVATE.")
    return ap.parse_args()


class Spider_ADC_MAP(ScanSpider):
    """Scan Spider: Spider_ADC_MAP

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param matlab_code: matlab code folder
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir,
                 xnat_project, xnat_subject, xnat_session, xnat_scan,
                 matlab_code,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix=""):
        """Entry point for Spider_ADC_MAP Class."""
        super(Spider_ADC_MAP,
              self).__init__(spider_path, jobdir, xnat_project, xnat_subject,
                             xnat_session, xnat_scan, xnat_host, xnat_user,
                             xnat_pass, suffix)
        self.inputs = list()
        self.matlab_code = matlab_code
        self.pdf_final = os.path.join(self.jobdir,
                                      '%s_ADC_MAPS_report.pdf'
                                      % xnat_session)

    def pre_run(self, argument_parse):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        resource = 'DICOM'  # resource to download from the scan on XNAT
        input_dir = XnatUtils.makedir(os.path.join(self.jobdir, 'inputs'))
        self.inputs.extend(self.download(self.xnat_scan, resource, input_dir))
        if len(self.inputs) == 1 and self.inputs[0].endswith('.zip'):
            self.time_writer('Unzipping the dicoms...')
            os.system('unzip -d %s -j %s > /dev/null'
                      % (os.path.join(input_dir, resource), self.inputs[0]))
            os.remove(self.inputs[0])
        self.inputs = get_dicom_list(input_dir)

    def run(self):
        """Method running the process for the spider on the inputs data."""
        output_dir = XnatUtils.makedir(os.path.join(self.jobdir, 'outputs'))
        osirix_dir = XnatUtils.makedir(os.path.join(self.jobdir, 'OsiriX'))
        self.time_writer('Dicom folder: %s' % os.path.dirname(self.inputs[0]))
        mat_lines = DEFAULT_ADC_TEMPLATE.format(
                matlab_code=self.matlab_code,
                input_path=os.path.dirname(self.inputs[0]),
                output=output_dir,
                pdf_name=self.pdf_final)
        matlab_script = os.path.join(output_dir, 'run_matlab_verdict.m')
        with open(matlab_script, "w") as f:
            f.writelines(mat_lines)
        self.run_matlab(matlab_script, verbose=True)

        # Zip outputs:
        # dcm_files = glob.glob(os.path.join(output_dir, '*', '*.dcm'))
        dcm_files = get_dicom_list(output_dir)
        for dicom in dcm_files:
            shutil.copy(dicom, osirix_dir)

        self.time_writer('Zipping OsiriX resource ...')
        # Zip the DICOMs output:
        initdir = os.getcwd()
        # Zip all the files in the directory
        zip_name = os.path.join(osirix_dir, 'osirix.zip')
        os.chdir(osirix_dir)
        os.system('zip -r %s * > /dev/null' % zip_name)
        # return to the initial directory:
        os.chdir(initdir)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        osirix_zip = os.path.join(self.jobdir, 'OsiriX', 'osirix.zip')
        results_dict = {'PDF': self.pdf_final,
                        'OsiriX': osirix_zip}

        self.upload_dict(results_dict)
        self.end()

    def run_matlab(self, matlab_script, verbose=False):
        """Call MATLAB with -nodesktop -nosplash and -singlecompthread.

        :param matlab_script: Full path to the .m file to run
        :param verbose: True to print all MATLAB output to terminal, False to
         suppress.
        :return: None
        """
        self.time_writer("Matlab script: %s running ..." % matlab_script)
        cmd = "matlab -nodisplay -nodesktop -nosplash -singleCompThread \
    < %s" % matlab_script
        if not verbose:
            matlabdir = os.path.dirname(matlab_script)
            prefix = os.path.basename(matlab_script).split('.')[0]
            cmd = cmd+' > '+os.path.join(matlabdir, prefix+'_outlog.log')
        os.system(cmd)
        self.time_writer("Matlab script: %s done" % matlab_script)


def is_dicom(fpath):
    """check if the file is a DICOM medical data.

    :param fpath: path of the file
    :return boolean: true if it's a DICOM, false otherwise
    """
    file_call = '''file {fpath}'''.format(fpath=fpath)
    output = sb.check_output(file_call.split())
    if 'dicom' in output.lower():
        return True

    return False


def get_dicom_list(directory):
    """get the list of DICOMs file from the directory and subdir.

    :param directory: directory containing the DICOM files.
    :return list(): list of filepaths that are dicoms in directory
    """
    fnames = os.listdir(directory)
    dicom_paths = list()
    for fname in fnames:
        fpath = os.path.join(directory, fname)
        if os.path.isfile(fpath) and is_dicom(fpath):
            dicom_paths.append(fpath)
        elif os.path.isdir(fpath):
            dicom_paths.extend(get_dicom_list(fpath))

    return dicom_paths


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_ADC_MAP(spider_path=sys.argv[0],
                                jobdir=args.temp_dir,
                                xnat_project=args.proj_label,
                                xnat_subject=args.subj_label,
                                xnat_session=args.sess_label,
                                xnat_scan=args.scan_label,
                                matlab_code=args.matlab_code,
                                xnat_host=args.host,
                                xnat_user=args.user,
                                xnat_pass=None,
                                suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault & Vasia Papoutsaki",
                          "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run(args)

    # Run method
    spider_obj.run()

    # Finish method to copy results
    if not args.skipfinish:
        spider_obj.finish()
