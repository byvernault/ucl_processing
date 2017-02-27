"""Spider_Stretched_Maps.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Stretched_Maps
Spider version: 1.0.0
Creation date:  2017-02-23 14:25:45
Purpose:        Calculate DDC maps and a maps using the DW Head and Neck data
"""

# Python packages import
import os
import sys

from dax import spiders, ScanSpider, XnatUtils


__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Calculate DDC maps and a maps using the DW Head and Neck data"
__spider_name__ = "Stretched_Maps"
__version__ = "1.0.0"
__modifications__ = """2017-02-23 14:25:45 - Original write"""


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
    ap = spiders.get_scan_argparser("Stretched_Maps", __purpose__)
    _h = "Matlab code folder containing XNAT_StretchedMaps_final.m"
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help=_h)
    return ap.parse_args()


class Spider_Stretched_Maps(ScanSpider):
    """Scan Spider: Spider_Stretched_Maps

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param matlab_code: matlab path to the code
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """
    def __init__(self, spider_path, jobdir,
                 xnat_project, xnat_subject, xnat_session, xnat_scan,
                 matlab_code,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix="", subdir=True, skip_finish=False):
        """Entry point for Spider_Stretched_Maps Class."""
        super(Spider_Stretched_Maps,
              self).__init__(spider_path, jobdir, xnat_project, xnat_subject,
                             xnat_session, xnat_scan, xnat_host, xnat_user,
                             xnat_pass, suffix, subdir, skip_finish)
        # Inputs to download from XNAT specified by:
        #   type: 'scan' or 'assessor' or 'session' or 'subject' or 'project'
        #   label: label on xnat for the object, e.g: '0002' for a scan
        #   resource: label of the resource on xnat, e.g: NIFTI
        #   dir (optional): directory where to download the data
        #   scan (optional): if using an scan assessor and giving just the
        #                    proctype to the label key, generate the
        #                    assessor_label string.
        self.inputs = [
            {'type': 'scan', 'label': xnat_scan, 'resource': 'DICOM'},
        ]
        _format = 'Report_Stretched_Maps_%s_%s.pdf'
        self.matlab_code = matlab_code
        self.pdf_final = os.path.join(
            self.jobdir, _format % (self.xnat_session, self.xnat_scan))

    def pre_run(self):
        """Method to download data from XNAT."""
        # Download inputs specified by self.inputs
        self.download_inputs()
        fname = os.path.basename(self.data[self.xnat_scan]['DICOM'][0])
        if fname == 'dicoms.zip':
            input_dir = os.path.dirname(self.data[self.xnat_scan]['DICOM'][0])
            XnatUtils.unzip_list(
                self.data[self.xnat_scan]['DICOM'][0], input_dir)

    def run(self):
        """Method running the process for the spider on the inputs data."""
        # Run command define by self.cmd_args
        matlab_template = """addpath('$mpath');
XNAT_StretchedMaps_final('$dcm_folder', '$ddc_maps', '$adc_maps', '$pdf');"""
        folder = os.path.dirname(self.data[self.xnat_scan]['DICOM'][0])
        dccdir = os.path.join(self.jobdir, 'dcc_maps')
        os.makedirs(dccdir)
        adcdir = os.path.join(self.jobdir, 'adc_maps')
        os.makedirs(adcdir)
        self.cmd_args = {
            'exe': 'matlab',
            'template': matlab_template,
            'args': {
                'dcm_folder': folder,
                'ddc_maps': dccdir,
                'adc_maps': adcdir,
                'mpath': self.matlab_code,
                'pdf': self.pdf_final,
            }
        }
        self.run_cmd_args()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        dccdir = os.path.join(self.jobdir, 'dcc_maps')
        adcdir = os.path.join(self.jobdir, 'adc_maps')
        results_dict = {'PDF': self.pdf_final,
                        'OsiriX': [adcdir, dccdir],
                        }
        self.upload_dict(results_dict)
        self.end()


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Stretched_Maps(
                    spider_path=sys.argv[0],
                    jobdir=args.temp_dir,
                    xnat_project=args.proj_label,
                    xnat_subject=args.subj_label,
                    xnat_session=args.sess_label,
                    xnat_scan=args.scan_label,
                    matlab_code=args.matlab_code,
                    xnat_host=args.host,
                    xnat_user=args.user,
                    xnat_pass=None,
                    suffix=args.suffix,
                    subdir=args.subdir,
                    skip_finish=args.skipfinish)

    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    if not args.skipfinish:
        spider_obj.finish()
