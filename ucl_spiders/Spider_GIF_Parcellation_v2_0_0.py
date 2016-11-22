"""Spider to perform gif parcellation.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    GIF_Parcellation
Spider version: 1.0.0
Creation date:  2016-02-05 13:52:08.737570
Purpose:        Parcellation of the brain using GIF: Geodesic Information Flow.
"""

# Python packages import
import os
import sys
import glob
from dax import spiders, ScanSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Parcellation of the brain using GIF: Geodesic Information Flow."
__spider_name__ = "GIF_Parcellation"
__version__ = "1.0.0"
__modifications__ = """2016-03-15 14:56 - Adding working_dir options
2016-05-10 18:04:01 - Update to new format respecting pep8
"""


GIF_CMD = """{exe_path} \
-i {input} \
-o {output} \
-d {db_xml} \
--no_qsub \
--openmp_core {number_core} \
--n_procs 1 \
--working_dir '{working_dir}' \
--remove_tmp
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
        --host   : host for XNAT (default: XNAT_HOST env variable)
        --user   : user on XNAT (default: XNAT_USER env variable)
    your arguments:
        --dbt    : gif-based database xml file describing the inputs
        --gif    : Path to the Gif python script: perform_gif_propagation.py
        --openmp_core : number of core use by reg_aladin. Default: one
        --working_dir : working directory for temp files
    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_scan_argparser("GIF_Parcellation", __purpose__)
    ap.add_argument("--dbt", dest="dbtemplate", required=True,
                    help="gif-based database xml file describing the inputs.")
    ap.add_argument("--gif", dest="gif_script", required=True,
                    help="Path to perform_gif_propagation.py.")
    ap.add_argument("--openmp_core", dest="openmp_core", default=1,
                    help="Number of core used by reg_aladin.")
    ap.add_argument("--working_dir", dest="working_dir", default="",
                    help="working directory for temp files. Default: output")
    return ap.parse_args()


class Spider_GIF_Parcellation(ScanSpider):
    """Spider for Parcellation of the brain using Geodesic Information Flow.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param number_core: number of core used by reg_aladin
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project,
                 xnat_subject, xnat_session, xnat_scan,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 number_core=1, suffix=""):
        """Entry point for Spider_GIF_Parcellation Class."""
        super(Spider_GIF_Parcellation, self).__init__(spider_path, jobdir,
                                                      xnat_project,
                                                      xnat_subject,
                                                      xnat_session,
                                                      xnat_scan,
                                                      xnat_host,
                                                      xnat_user,
                                                      xnat_pass,
                                                      suffix, subdir=False)
        self.inputs = list()
        self.number_core = number_core
        # Print version for Niftyreg - GIFi
        self.check_executable('reg_aladin', 'reg_aladin')
        self.pdf_final = os.path.join(self.jobdir, 'GIF_parcellation.pdf')

    def pre_run(self):
        """Method to download data from XNAT.

        :return: None
        """
        resource = 'NIFTI'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        if not os.path.exists(os.path.join(self.jobdir, 'inputs')):
            os.makedirs(folder)
        self.inputs.extend(self.download(self.xnat_scan, resource, folder))

    def run(self, gif_script, dbtemplate, working_dir):
        """Method running the process for the spider on the inputs data.

        :param gif_script: Path to the Gif python script
        :param dbtemplate: gif-based database xml file describing the inputs
        :return: None
        """
        if len(working_dir) > 0 and not os.path.exists(working_dir):
            os.makedirs(working_dir)

        if gif_script.endswith('perform_gif_propagation.py'):
            exe_path = gif_script
        else:
            exe_path = os.path.join(gif_script, "perform_gif_propagation.py")

        if not os.path.exists(exe_path):
            raise Exception("Python script: %s not found" % (exe_path))
        else:
            if not os.path.exists(os.path.join(self.jobdir, 'outputs')):
                os.makedirs(os.path.join(self.jobdir, 'outputs'))
            cmd = GIF_CMD.format(exe_path=exe_path,
                                 input=self.inputs[0],
                                 output=os.path.join(self.jobdir, 'outputs'),
                                 db_xml=dbtemplate,
                                 number_core=self.number_core,
                                 working_dir=working_dir)
            self.run_system_cmd(cmd)
            self.make_pdf()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        out_dir = os.path.join(self.jobdir, 'outputs')
        # Images outputs:
        bias_corrected = glob.glob(os.path.join(out_dir,
                                                '*bias_corrected.nii.gz'))
        brain = glob.glob(os.path.join(out_dir, '*brain.nii.gz'))
        labels = glob.glob(os.path.join(out_dir, '*labels.nii.gz'))
        prior = glob.glob(os.path.join(out_dir, '*prior.nii.gz'))
        seg = glob.glob(os.path.join(out_dir, '*seg.nii.gz'))
        tiv = glob.glob(os.path.join(out_dir, '*tiv.nii.gz'))
        # Volumes:
        volumes = glob.glob(os.path.join(out_dir, '*volumes.xml'))
        results_dict = {'PDF': self.pdf_final,
                        'BIAS_COR': bias_corrected,
                        'BRAIN': brain,
                        'LABELS': labels,
                        'PRIOR': prior,
                        'SEG': seg,
                        'TIV': tiv,
                        'STATS': volumes}
        self.upload_dict(results_dict)
        self.end()

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        # PDF pages:
        pdf_pages = {
            '1': os.path.join(self.jobdir, 'GIF_parcellation_page1.pdf'),
            '2': os.path.join(self.jobdir, 'GIF_parcellation_page2.pdf')
        }

        # Images outputs:
        out_dir = os.path.join(self.jobdir, 'outputs')
        bias_corrected = glob.glob(os.path.join(out_dir,
                                                '*bias_corrected.nii.gz'))
        brain = glob.glob(os.path.join(out_dir, '*brain.nii.gz'))
        labels = glob.glob(os.path.join(out_dir, '*labels.nii.gz'))
        prior = glob.glob(os.path.join(out_dir, '*prior.nii.gz'))
        seg = glob.glob(os.path.join(out_dir, '*seg.nii.gz'))
        tiv = glob.glob(os.path.join(out_dir, '*tiv.nii.gz'))
        list_images = [bias_corrected, brain, labels, seg, tiv, prior]

        # Page 1:
        images = []
        for index, image_file in enumerate(list_images):
            if len(image_file) != 1:
                err = '%s output image not found or more than one file found.'
                raise Exception(err % (image_file))
            images.append(image_file[0])

        labels = {
            '0': 'Bias Corrected',
            '1': 'Brain',
            '2': 'Labels',
            '3': 'Segmentation',
            '4': 'tiv',
            '5': 'prior'
        }
        cmap = {
            '0': 'gray',
            '1': 'gray',
            '2': None,
            '3': 'gray',
            '4': 'gray',
            '5': None
        }
        self.plot_images_page(pdf_pages['1'], 1, images,
                              'GIF_Parcellation Pipeline',
                              image_labels=labels, cmap=cmap)

        # Page 2
        """ # Volumes:
        volumes = glob.glob(os.path.join(out_dir, '*volumes.xml'))
        if len(volumes) != 1:
            err = '%s output csv file with information on volumes not found \
or more than one file found.'
            raise Exception(err % (volumes))
        with open(volumes[0], 'rb') as csvfileread:
            csvreader = csv.reader(csvfileread, delimiter=',')
            li_name = csvreader.next()
            li_volume = csvreader.next()

        di_stats = OrderedDict()
        for index, name in enumerate(li_name):
            di_stats[name] = li_volume[index]

        self.plot_stats_page(pdf_pages['2'], 2, di_stats,
                             'Volumes computed by GIF_Parcellation',
                             columns_header=['Label Name', 'Volume'])

        # Join the two pages for the PDF:
        self.merge_pdf_pages(pdf_pages, self.pdf_final)
        """
        os.rename(pdf_pages['1'], self.pdf_final)


if __name__ == '__main__':
    # arguments
    args = parse_args()

    # generate spider object:
    spider_obj = Spider_GIF_Parcellation(spider_path=sys.argv[0],
                                         jobdir=args.temp_dir,
                                         xnat_project=args.proj_label,
                                         xnat_subject=args.subj_label,
                                         xnat_session=args.sess_label,
                                         xnat_scan=args.scan_label,
                                         xnat_host=args.host,
                                         xnat_user=args.user,
                                         xnat_pass=None,
                                         number_core=args.openmp_core,
                                         suffix=args.suffix)

    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run(args.gif_script, args.dbtemplate, args.working_dir)

    # Finish method to copy results
    spider_obj.finish()
