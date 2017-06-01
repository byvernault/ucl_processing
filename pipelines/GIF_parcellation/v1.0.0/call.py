import os
import csv
import glob
from dax import spiders, XnatUtils
from collections import OrderedDict


JOB_DIR = '${temp_dir}'
IN_FILE = '${t1_file}'
DB_TEMPLATE = '${dbt}'
GIF_EXE = '${gif}'
OPENMP_CORE = '${openmp_core}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''{exe_path} \
-i {input} \
-o {output} \
-d {db_xml} \
--no_qsub \
--openmp_core {number_core} \
--n_procs 1 \
--working_dir '{working_dir}' \
--remove_tmp'''


def main():
    if len(WORKING_DIR) > 0 and not os.path.exists(WORKING_DIR):
        os.makedirs(WORKING_DIR)

    if os.path.exists(GIF_EXE) or XnatUtils.executable_exists(GIF_EXE):
        cmd = EXE_CMD.format(exe_path=GIF_EXE,
                             input=IN_FILE,
                             output=JOB_DIR,
                             db_xml=DB_TEMPLATE,
                             number_core=OPENMP_CORE,
                             working_dir=WORKING_DIR)
        os.system(cmd)
        make_pdf()
    else:
        raise Exception("Error: %s not found" % (GIF_EXE))


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
    # Volumes:
    volumes = glob.glob(os.path.join(JOB_DIR, '*volumes.csv'))

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
    spiders.plot_images(pdf_pages['1'], 1, images,
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

    spiders.plot_stats(pdf_pages['2'], 2, di_stats,
                       'Volumes computed by GIF_Parcellation',
                       columns_header=['Label Name', 'Volume'])

    # Join the two pages for the PDF:
    pdf_final = os.path.join(JOB_DIR, 'GIF_parcellation.pdf')
    spiders.merge_pdf(pdf_pages, pdf_final)


if __name__ == '__main__':
    main()
