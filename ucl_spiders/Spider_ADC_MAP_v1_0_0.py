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
from dax import XnatUtils, spiders, ScanSpider, SessionSpider

__author__ = "Benjamin Yvernault & Vasia Papoutsaki"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Compute ADC Maps for HandN project"
__spider_name__ = "ADC_MAP"
__version__ = "1.0.0"
__modifications__ = """2016-12-01 13:49:23.748867 - Original write"""


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
        resource = 'XXXX'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        self.inputs.extend(self.download(self.xnat_scan, resource, folder))

        #
        # YOU CAN ADD MORE line to download other data LIKE
        #    self.inputs.extend(self.download(**kwargs))
        #

    def run(self):
        """Method running the process for the spider on the inputs data."""
        #
        # CODE THAT YOU WANT TO RUN ON THE DATA
        #

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': 'pdfpath'
                        #
                        # ADD OTHER RESULTS YOU WANT TO SAVE
                        #
                        }
        self.upload_dict(results_dict)
        self.end()


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_ADC_MAP(spider_path=sys.argv[0],
                               jobdir=args.temp_dir,
                               xnat_project=args.proj_label,
                               xnat_subject=args.subj_label,
                               xnat_session=args.sess_label,
                               xnat_scan=args.scan_label,
                               xnat_host=args.host,
                               xnat_user=args.user,
                               xnat_pass=None,
                               suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault & Vasia Papoutsaki", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run(args)

    # Run method
    spider_obj.run()

    # Finish method to copy results
    spider_obj.finish()
