import os
import glob
import xml.etree.ElementTree as ET
from dax import spiders, XnatUtils
from collections import OrderedDict


JOB_DIR = '${temp_dir}'
IN_FILE = '${t1_file}'
DB_TEMPLATE = '${dbt}'
NIFTYPIPE_EXE = '${gif}'
OPENMP_CORE = '${openmp_core}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''{exe} \
-i {input} \
-o {output} \
-d {db_xml} \
--no_qsub \
--n_procs 1 \
--remove_tmp \
{omp} \
{wdir}'''
OMP = '''--openmp_core {number_core}'''
WDIR = '''--working_dir '{working_dir}' '''


def main():
    """ Main function."""
    _working_dir = None
    if WORKING_DIR != 'None':
        _working_dir = os.path.join(WORKING_DIR, '${assessor_label}')
        if not os.path.exists(_working_dir):
            os.makedirs(_working_dir)

    if os.path.exists(NIFTYPIPE_EXE) or \
       XnatUtils.executable_exists(NIFTYPIPE_EXE):
        if OPENMP_CORE is not None and OPENMP_CORE != 'None':
            _omp = OMP.format(number_core=OPENMP_CORE)
        if _working_dir is not None:
            _wd = WDIR.format(working_dir=_working_dir)
        cmd = EXE_CMD.format(exe=NIFTYPIPE_EXE,
                             input=IN_FILE,
                             output=JOB_DIR,
                             db_xml=DB_TEMPLATE,
                             omp=_omp, wdir=_wd)
        os.system(cmd)
        make_pdf()
    else:
        raise Exception("Error: %s not found" % (NIFTYPIPE_EXE))


def make_pdf():
    """Method to make the PDF for the spider.

    :return: None
    """
    # PDF pages:
    pdf_pages = {
        '1': os.path.join(JOB_DIR, 'GIF_parcellation_page1.pdf'),
        '2': os.path.join(JOB_DIR, 'GIF_parcellation_page2.pdf')
    }

    # Images outputs:
    bias_corrected = glob.glob(os.path.join(JOB_DIR,
                                            '*bias_corrected.nii.gz'))
    brain = glob.glob(os.path.join(JOB_DIR, '*brain.nii.gz'))
    labels = glob.glob(os.path.join(JOB_DIR, '*labels.nii.gz'))
    prior = glob.glob(os.path.join(JOB_DIR, '*prior.nii.gz'))
    seg = glob.glob(os.path.join(JOB_DIR, '*seg.nii.gz'))
    tiv = glob.glob(os.path.join(JOB_DIR, '*tiv.nii.gz'))
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
    spiders.plot_images_page(pdf_pages['1'], 1, images,
                             'GIF_Parcellation Pipeline',
                             image_labels=labels, cmap=cmap)

    # Page 2
    # Volumes:
    volumes = glob.glob(os.path.join(JOB_DIR, '*volumes.xml'))
    if len(volumes) != 1:
        err = '%s output csv file with information on volumes not found \
or more than one file found.'
        raise Exception(err % (volumes))
    tree = ET.parse(volumes[0])
    root = tree.getroot()
    di_stats = OrderedDict()
    for tissue in root.findall('tissues'):
        for item in tissue.findall('item'):
            di_stats[item.find('name').text] = item.find('volumeProb').text
    for tissue in root.findall('labels'):
        for item in tissue.findall('item'):
            di_stats[item.find('name').text] = item.find('volumeProb').text

    spiders.plot_stats(pdf_pages['2'], 2, di_stats,
                       'Volumes computed by GIF_Parcellation',
                       columns_header=['Label Name', 'Volume ml'])

    # Join the two pages for the PDF:
    pdf_final = os.path.join(JOB_DIR, 'GIF_parcellation.pdf')
    spiders.merge_pdfs(pdf_pages, pdf_final)


if __name__ == '__main__':
    main()
