from dax import XnatUtils
import glob
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.misc import imread


JOB_DIR = '${temp_dir}'
ENV_SOURCE = '${env_source}'
DWIS = '${dwi}'
BVALS = '${bval}'
BVECS = '${bvec}'
NIFTYPIPE_EXE = '${exe}'
OPENMP_CORE = '${omp}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''${exe} \
-i {dwis} \
-a {bvals} \
-e {bvecs} \
-t ${t1} \
--t1_mask ${tiv} \
-o ${temp_dir} \
--no_qsub \
--n_procs 1 \
--rot 34.56 \
--etd 2.46 \
--ped=-y \
--remove_tmp \
{omp} \
{wdir}'''
OMP = '''--openmp_core {number_core}'''
WDIR = '''--working_dir '{working_dir}' '''


def main():
    """ Main function."""
    if ENV_SOURCE is not None and os.path.isfile(ENV_SOURCE):
        os.system('sh {}'.format(ENV_SOURCE))

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

        cmd = EXE_CMD.format(
            dwis=' '.join(DWIS.split(',')),
            bvals=' '.join(BVALS.split(',')),
            bvecs=' '.join(BVECS.split(',')),
            omp=_omp, wdir=_wd)
        os.system(cmd)
        make_pdf()
    else:
        raise Exception("Error: executable %s not found" % (NIFTYPIPE_EXE))


def make_pdf():
    """Method to make the PDF for the spider.

    :return: None
    """
    # PDF pages:
    pdf_report = os.path.join(JOB_DIR, 'DWI_2_DTI_report.pdf')

    # Images outputs:
    png_qc = glob.glob(os.path.join(JOB_DIR, '*.png'))

    if len(png_qc) < 2:
        raise Exception('Missing PNG files for PDF.')

    pp = PdfPages(pdf_report)
    plt.figure(1)
    plot_image(png_qc[0])
    pp.savefig(plt.gcf())  # This generates page 1
    plt.figure(2)
    plot_image(png_qc[1])
    pp.savefig(plt.gcf())  # This generates page 2
    pp.close()


def plot_image(filepath):
    im = imread(filepath).astype(np.float32) / 255
    plt.imshow(im)
    a = plt.gca()
    a.get_xaxis().set_visible(False)  # We don't need axis ticks
    a.get_yaxis().set_visible(False)


if __name__ == '__main__':
    main()
