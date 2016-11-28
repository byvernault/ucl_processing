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
import csv
import glob
import numpy as np
import nibabel as nib
from datetime import datetime
import matplotlib.pyplot as plt
from collections import OrderedDict
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
--remove_tmp"""


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
        volumes = glob.glob(os.path.join(out_dir, '*volumes.csv'))
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
        # Volumes:
        volumes = glob.glob(os.path.join(out_dir, '*volumes.csv'))

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

    def plot_images_page(self, pdf_path, page_index, nii_images, title,
                         image_labels, slices=None, cmap='gray',
                         vmins=None, vmaxs=None):
        """Plot list of images (3D-4D) on a figure (PDF page).
        plot_images_figure will create one pdf page with only images.
        Each image corresponds to one line with by default axial/sag/cor view
        of the mid slice. If you use slices, it will show different slices of
        the axial plan view. You can specify the cmap and the vmins/vmaxs if
        needed by using a dictionary with the index of each line (0, 1, ...).
        :param pdf_path: path to the pdf to save this figure to
        :param page_index: page index for PDF
        :param nii_images: python list of nifty images
        :param title: Title for the report page
        :param image_labels: list of titles for each images
            one per image in nii_images
        :param slices: dictionary of list of slices to display
            if None, display axial, coronal, sagital
        :param cmap: cmap to use to display images or dict
            of cmaps for each images with the indices as key
        :param vmins: define vmin for display (dict)
        :param vmaxs: define vmax for display (dict)
        :return: pdf path created
        E.g for two images:
        images = [imag1, image2]
        slices = {'0':[50, 80, 100, 130],
                  '1':[150, 180, 200, 220]}
        labels = {'0': 'Label 1',
                  '1': 'Label 2'}
        cmaps = {'0':'hot',
                 '1': 'gray'}
        vmins = {'0':10,
                 '1':20}
        vmaxs = {'0':100,
                 '1':150}
        """
        plt.ioff()
        self.time_writer('INFO: generating pdf page %d with images.'
                         % page_index)
        fig = plt.figure(page_index, figsize=(7.5, 10))
        # Titles:
        if not isinstance(cmap, dict):
            default_cmap = cmap
            cmap = {}
        else:
            default_cmap = 'gray'
        if not isinstance(vmins, dict):
            self.time_writer("Warning: vmins wasn't a dictionary. \
Using default.")
            vmins = {}
        if not isinstance(vmaxs, dict):
            self.time_writer("Warning: vmaxs wasnt' a dictionary. \
Using default.")
            vmaxs = {}
        if isinstance(nii_images, str):
            nii_images = [nii_images]
        number_im = len(nii_images)
        for index, image in enumerate(nii_images):
            # Open niftis with nibabel
            f_img = nib.load(image)
            data = f_img.get_data()
            if len(data.shape) == 4:
                data = data[:, :, :, data.shape[3]/2]
            default_slices = [data.shape[2]/4, data.shape[2]/2,
                              3*data.shape[2]/4]
            default_label = 'Line %s' % index

            if slices:
                if not isinstance(slices, dict):
                    self.time_writer("Warning: slices wasn't a dictionary. \
Using default.")
                    slices = {}
                self.time_writer('INFO: showing different slices.')
                li_slices = slices.get(str(index), default_slices)
                slices_number = len(li_slices)
                for slice_ind, slice_value in enumerate(li_slices):
                    ind = slices_number*index+slice_ind+1
                    ax = fig.add_subplot(number_im, slices_number, ind)
                    data_z_rot = np.rot90(data[:, :, slice_value])
                    ax.imshow(data_z_rot,
                              cmap=cmap.get(str(index), default_cmap),
                              vmin=vmins.get(str(index), None),
                              vmax=vmaxs.get(str(index), None))
                    ax.set_title('Slice %d' % slice_value, fontsize=7)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    if slice_ind == 0:
                        ax.set_ylabel(image_labels.get(str(index),
                                      default_label), fontsize=9)
            else:
                self.time_writer('INFO: display different plan view \
(ax/sag/cor) of the mid slice.')
                ax = fig.add_subplot(number_im, 3, 3*index+1)
                data_z_rot = np.rot90(data[:, :, data.shape[2]/2])
                ax.imshow(data_z_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Axial', fontsize=7)
                ax.set_ylabel(image_labels.get(str(index), default_label),
                              fontsize=9)
                ax.set_xticks([])
                ax.set_yticks([])
                ax = fig.add_subplot(number_im, 3, 3*index+2)
                data_y_rot = np.rot90(data[:, data.shape[1]/2, :])
                ax.imshow(data_y_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Coronal', fontsize=7)
                ax.set_axis_off()
                ax = fig.add_subplot(number_im, 3, 3*index+3)
                data_x_rot = np.rot90(data[data.shape[0]/2, :, :])
                ax.imshow(data_x_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Sagittal', fontsize=7)
                ax.set_axis_off()

        fig.tight_layout()
        date = datetime.now()
        # Titles page
        plt.figtext(0.5, 0.985, '-- %s PDF report --' % title,
                    horizontalalignment='center', fontsize=12)
        plt.figtext(0.5, 0.02, 'Date: %s -- page %d' % (str(date), page_index),
                    horizontalalignment='center', fontsize=8)
        fig.savefig(pdf_path, transparent=True, orientation='portrait',
                    dpi=100)
        plt.close(fig)
        return pdf_path

    def plot_stats_page(self, pdf_path, page_index, stats_dict, title,
                        tables_number=3, columns_header=['Header', 'Value'],
                        limit_size_text_column1=30,
                        limit_size_text_column2=10):
        """Generate pdf report of stats information from a csv/txt.
        plot_stats_page generate a pdf page displaying a dictionary
        of stats given to the function. Column 1 represents the key
        or header and the column 2 represents the value associated.
        You can rename the two column by using the args column1/2.
        There are three columns than can have 50 values max.
        :param pdf_path: path to the pdf to save this figure to
        :param page_index: page index for PDF
        :param stats_dict: python dictionary of key=value to display
        :param title: Title for the report page
        :param tables_number: number of columns to display (def:3)
        :param columns_header: list of header for the column
            default: header, value
        :param limit_size_text_column1: limit of text display in column 1
        :param limit_size_text_column2: limit of text display in column 2
        :return: pdf path created
        """
        plt.ioff()
        self.time_writer('INFO: generating pdf page %d with stats.'
                         % page_index)
        cell_text = list()
        for key, value in stats_dict.items():
            txt = smaller_str(key.strip().replace('"', ''),
                              size=limit_size_text_column1)
            val = smaller_str(str(value),
                              size=limit_size_text_column2)
            cell_text.append([txt, "%s" % val])

        # Make the table
        fig = plt.figure(page_index, figsize=(7.5, 10))
        nb_stats = len(stats_dict.keys())
        for i in range(tables_number):
            ax = fig.add_subplot(1, tables_number, i+1)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.axis('off')
            the_table = ax.table(
                    cellText=cell_text[nb_stats/3*i:nb_stats/3*(i+1)],
                    colColours=[(0.8, 0.4, 0.4), (1.0, 1.0, 0.4)],
                    colLabels=columns_header,
                    colWidths=[0.8, 0.32],
                    loc='center',
                    rowLoc='left',
                    colLoc='left',
                    cellLoc='left')

            the_table.auto_set_font_size(False)
            the_table.set_fontsize(6)

        # Set footer and title
        date = datetime.now()
        plt.figtext(0.5, 0.985, '-- %s PDF report --' % title,
                    horizontalalignment='center', fontsize=12)
        plt.figtext(0.5, 0.02, 'Date: %s -- page %d' % (str(date), page_index),
                    horizontalalignment='center', fontsize=8)
        fig.savefig(pdf_path, transparent=True, orientation='portrait',
                    dpi=300)
        plt.close(fig)
        return pdf_path

    def merge_pdf_pages(self, pdf_pages, pdf_final):
        """Concatenate all pdf pages in the list into a final pdf.
        You can provide a list of pdf path or give a dictionary
        with each page specify by a number:
          pdf_pages = {'1': pdf_page1, '2': pdf_page2}
        :param pdf_pages: python list or dictionary of pdf page path
        :param pdf_final: final PDF path
        :return: pdf path created
        """
        self.time_writer('INFO: Concatenate all pdfs pages.')
        pages = ''
        if isinstance(pdf_pages, dict):
            for key in sorted(pdf_pages.iterkeys()):
                pages += '%s %s ' % (pages, pdf_pages[key])
        elif isinstance(pdf_pages, list):
            pages = ' '.join(pdf_pages)
        else:
            raise Exception('Wrong type for pdf_pages (list or dict).')
        cmd = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH \
-sDEVICE=pdfwrite -sOutputFile=%s %s' % (pdf_final, pages)
        self.time_writer('INFO:saving final PDF: %s ' % cmd)
        self.run_system_cmd(cmd)
        return pdf_final


def smaller_str(str_option, size=10, end=False):
    """Method to shorten a string into a smaller size.
    :param str_option: string to shorten
    :param size: size of the string to keep (default: 10 characters)
    :param end: keep the end of the string visible (default beginning)
    :return: shortened string
    """
    if len(str_option) > size:
        if end:
            return '...%s' % (str_option[-size:])
        else:
            return '%s...' % (str_option[:size])
    else:
        return str_option


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
