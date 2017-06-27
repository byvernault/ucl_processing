import os
import glob
from dax import spiders, XnatUtils


JOB_DIR = '${temp_dir}'
ASSESSOR = '${assessor_label}'
GAD_FILE = '${gad}'
MPRAGE_FILE = '${mprage}'
LABELS_FILE = '${labels}'
BRAIN_FILE = '${brain}'
NIFTYPIPE_EXE = '${exe}'
OPENMP_CORE = '${omp}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''{exe_path} \
--ref {gad} \
--in_file {mprage} \
--labels {labels} \
--brain {brain} \
--output {output} \
--no_qsub \
--remove_tmp \
--n_procs 1 \
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

    if XnatUtils.executable_exists(NIFTYPIPE_EXE):
        _omp = ''
        _wd = ''
        if OPENMP_CORE is not None and OPENMP_CORE != 'None':
            _omp = OMP.format(number_core=OPENMP_CORE)
        if _working_dir is not None:
            _wd = WDIR.format(working_dir=_working_dir)

        cmd = EXE_CMD.format(exe_path=NIFTYPIPE_EXE,
                             gad=GAD_FILE,
                             mprage=MPRAGE_FILE,
                             labels=LABELS_FILE,
                             brain=BRAIN_FILE,
                             output=JOB_DIR,
                             omp=_omp,
                             wdir=_wd)
        os.system(cmd)
        rename_outputs()
        make_pdf()
    else:
        raise Exception("Error: executable %s not found" % (NIFTYPIPE_EXE))


def make_pdf():
    """Function to make the PDF for the spider.

    :return: None
    """
    # PDF pages:
    pdf = os.path.join(JOB_DIR, 'mprage2gad.pdf')

    # Images outputs:
    images = [GAD_FILE]
    images.extend(glob.glob(os.path.join(JOB_DIR, '*.nii.gz')))

    labels = dict()
    for ind in range(len(images)):
        if ind == 0:
            labels[str(ind)] = 'Reference'
        else:
            labels[str(ind)] = os.path.basename(
                images[ind]).split('2gad_registered')[0]

    spiders.plot_images(pdf, 1, images, 'MPRAGE to Gad Registration Pipeline',
                        image_labels=labels)


def rename_outputs():
    """Function to rename the outputs.

    :return: None
    """
    assr_handler = XnatUtils.AssessorHandler(ASSESSOR)
    session = assr_handler.get_session_label()
    _ftmp = '{}_{}2gad_registered.nii.gz'
    for image in glob.glob(os.path.join(JOB_DIR, '*.nii.gz')):
        filename = os.path.basename(image)
        if 'brain' in filename:
            fname = os.path.join(JOB_DIR, _ftmp.format(session, 'brain'))
        elif 'labels' in filename:
            fname = os.path.join(JOB_DIR, _ftmp.format(session, 'labels'))
        else:
            fname = os.path.join(JOB_DIR, _ftmp.format(session, 'mprage'))
        os.rename(image, fname)


if __name__ == '__main__':
    main()