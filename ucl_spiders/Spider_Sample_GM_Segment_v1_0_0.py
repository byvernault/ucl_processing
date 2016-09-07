"""Spider_Sample_GM_Segment.

Author:         Benjamin Yvernault & Dave Cash
contact:        b.yvernault@ucl.ac.uk
Spider name:    Sample_GM_Segment
Spider version: 1.0.0
Creation date:  2016-09-02 15:59:57.932632
Purpose:        Run sample GM segment .m file using SPM12
"""

# Python packages import
import os
import sys
import nibabel as nib
from dax import XnatUtils, spiders, ScanSpider

__author__ = "Benjamin Yvernault & Dave Cash"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Run sample GM segment .m file using SPM12"
__spider_name__ = "Sample_GM_Segment"
__version__ = "1.0.0"
__modifications__ = """2016-09-02 15:59:57.932632 - Original write"""

MAT_TEMPLATE = """
addpath('{matlab_code}');
sample_GM_segment('{input_file}', '{spm12}');
"""


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
    ap = spiders.get_scan_argparser("Sample_GM_Segment", __purpose__)
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help="Matlab code path to run sample GM segment.")
    ap.add_argument("--spm12", dest="spm12", default=None, required=True,
                    help="SPM 12 path.")
    return ap.parse_args()


class Spider_Sample_GM_Segment(ScanSpider):
    """Scan Spider: Spider_Sample_GM_Segment

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param mc: matlab code path
    :param spm12: SPM 12 path
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir,
                 xnat_project, xnat_subject, xnat_session, xnat_scan,
                 matlab_code, spm12,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix=""):
        """Entry point for Spider_Sample_GM_Segment Class."""
        super(Spider_Sample_GM_Segment,
              self).__init__(spider_path, jobdir, xnat_project, xnat_subject,
                             xnat_session, xnat_scan, xnat_host, xnat_user,
                             xnat_pass, suffix)
        self.inputs = list()
        self.input_file = ''
        self.matlab_code = os.path.abspath(matlab_code)
        self.spm12 = os.path.abspath(spm12)
        self.pdf_final = os.path.join(self.jobdir,
                                      '%s_sample_gm_report.pdf' % xnat_session)

    def pre_run(self, argument_parse):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        resource = 'NIFTI'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'Sample_GM_Segment')
        if not os.path.exists(os.path.join(self.jobdir, 'inputs')):
            os.makedirs(folder)
        self.inputs.extend(self.download(self.xnat_scan, resource, folder))

    def run(self):
        """Method running the process for the spider on the inputs data."""
        # Gzip files
        input_file = ''
        for filepath in self.inputs:
            if filepath.endswith('.nii.gz'):
                XnatUtils.gunzip_file(filepath)
                input_file = filepath[:-3]
        self.input_file = input_file
        folder = os.path.join(self.jobdir, 'Sample_GM_Segment')
        mat_lines = MAT_TEMPLATE.format(
                        matlab_code=self.matlab_code,
                        input_file=input_file,
                        spm12=self.spm12)
        matlab_script = os.path.join(folder, 'run_sample_GM.m')
        with open(matlab_script, "w") as f:
            f.writelines(mat_lines)
        XnatUtils.run_matlab(matlab_script, verbose=True)

        # Make report:
        self.make_pdf()

        # Gzip nii:
        XnatUtils.gzip_nii(folder)

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        folder = os.path.dirname(self.input_file)
        in_filename = os.path.basename(self.input_file)
        pdf_title = 'Sample GM Segment - SPM12 - Report'
        f_img = nib.load(os.path.join(folder, 'c1%s' % in_filename))
        data = f_img.get_data()
        nb_slices = data.shape[2]
        slices = {'0': [(nb_slices-1)/4, 2*(nb_slices-1)/5,
                        3*(nb_slices-1)/5, 3*(nb_slices-1)/4],
                  '1': [(nb_slices-1)/4, 2*(nb_slices-1)/5,
                        3*(nb_slices-1)/5, 3*(nb_slices-1)/4],
                  '2': [(nb_slices-1)/4, 2*(nb_slices-1)/5,
                        3*(nb_slices-1)/5, 3*(nb_slices-1)/4],
                  '3': [(nb_slices-1)/4, 2*(nb_slices-1)/5,
                        3*(nb_slices-1)/5, 3*(nb_slices-1)/4]}
        labels = {'0': 'c1',
                  '1': 'c2',
                  '2': 'c3',
                  '3': 'm'}
        images = [os.path.join(folder, 'c1%s' % in_filename),
                  os.path.join(folder, 'c2%s' % in_filename),
                  os.path.join(folder, 'c3%s' % in_filename),
                  os.path.join(folder, 'm%s' % in_filename)]
        self.plot_images_page(self.pdf_final, 1, images,
                              pdf_title, slices=slices, image_labels=labels)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        # Groupe Registrations:
        folder = os.path.dirname(self.input_file)
        rc_res = list()
        c_res = list()
        bias_res = list()
        mat_res = list()
        for filename in os.listdir(folder):
            if filename.startswith('c'):
                c_res.append(os.path.join(folder, filename))
            if filename.startswith('m'):
                bias_res.append(os.path.join(folder, filename))
            if filename.endswith('.mat'):
                mat_res.append(os.path.join(folder, filename))
            if filename.startswith('rc'):
                rc_res.append(os.path.join(folder, filename))
        mat_res.append(os.path.join(self.jobdir, 'Sample_GM_Segment',
                                    'run_sample_GM.m'))
        results_dict = {'PDF': self.pdf_final,
                        'RC': rc_res,
                        'BIAS': bias_res,
                        'C': c_res,
                        'MAT': mat_res
                        }
        self.upload_dict(results_dict)
        self.end()


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Sample_GM_Segment(
                               spider_path=sys.argv[0],
                               jobdir=args.temp_dir,
                               xnat_project=args.proj_label,
                               xnat_subject=args.subj_label,
                               xnat_session=args.sess_label,
                               xnat_scan=args.scan_label,
                               matlab_code=args.matlab_code,
                               spm12=args.spm12,
                               xnat_host=args.host,
                               xnat_user=args.user,
                               xnat_pass=None,
                               suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault & Dave Cash",
                          "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run(args)

    # Run method
    spider_obj.run()

    # Finish method to copy results
    spider_obj.finish()
