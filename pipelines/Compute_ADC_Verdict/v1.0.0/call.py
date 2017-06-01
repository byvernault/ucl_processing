import datetime
from dax import XnatUtils, spiders
import dicom
from dicom.dataset import Dataset, FileDataset
import glob
import numpy as np
import nibabel as nib
import os
import shutil


# Inputs:
LABEL = '${assessor_label}'
JOBDIR = '${temp_dir}'
ACQ_FILES = '${acquisitions}'.split(',')
DCM_FILE = '${dcm_file}'.split(',')[0]

# Variables
DEFAULT_MODEL = "VerdictProstate_Rmaps"
DEFAULT_VERDICT_TEMPLATE = """
addpath(genpath('{matlab_code}'));
compute_ADC_VERDICT('{input_path}',\
'{output}',\
'{scheme_filename}',\
'{camino}');
"""


def main():
    """ Main function."""
    if len(JOBDIR) > 0 and not os.path.exists(JOBDIR):
        os.makedirs(JOBDIR)

    dcm_folder = XnatUtils.makedir(os.path.join(JOBDIR, 'OsiriX'))

    project, subject, session = LABEL.split('-x-')[:3]
    for nb_acq, acq_path in enumerate(ACQ_FILES):
        nb_acq = nb_acq + 1
        if not os.path.exists(acq_path):
            continue
        else:
            folder = os.path.join(JOBDIR, str(nb_acq))
            os.makedirs(folder)

            # Unzip niftis:
            XnatUtils.gunzip_file(acq_path)
            acq_path = acq_path[:-3]

            # Command
            mat_lines = DEFAULT_VERDICT_TEMPLATE.format(
                matlab_code='${matlab_src}',
                input_path=acq_path,
                output=folder,
                scheme_filename='${scheme_filename}',
                camino='${camino}'
            )
            matlab_script = os.path.join(JOBDIR,
                                         'run_adc_map_verdict%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            cmd = "matlab -nodisplay -nodesktop -nosplash -singleCompThread < "
            os.system(cmd + matlab_script)

            # Generate Dicom for OsiriX
            res_nii = os.path.join(folder, 'FIT_ADC.nii')
            out_nii = os.path.join(folder,
                                   '%s_FIT_ADC_%d.nii' % (session, nb_acq))
            shutil.move(res_nii, out_nii)

            # Load dicom headers
            if not os.path.isfile(DCM_FILE):
                err = "DICOM File %s not found."
                raise Exception(err % DCM_FILE)
            sour_obj = dicom.read_file(DCM_FILE)

            # Convert all niftis to dicoms
            convert_nifti_2_dicoms(
                out_nii,
                sour_obj,
                dcm_folder,
                nb_acq,
                project,
                subject,
                session
            )

            # Gzip nii:
            XnatUtils.gzip_nii(folder)

    # Make pdf:
    make_pdf(session)

    # Remove reorient images:
    reorient_images = glob.glob(os.path.join(JOBDIR, '*', '*_reorient.nii.gz'))
    for reorient_im in reorient_images:
        os.remove(reorient_im)

    # Zip the DICOMs output:
    initdir = os.getcwd()
    # Zip all the files in the directory
    zip_name = os.path.join(JOBDIR, 'OsiriX', 'osirix.zip')
    os.chdir(os.path.join(JOBDIR, 'OsiriX'))
    os.system('zip -r %s * > /dev/null' % zip_name)

    # return to the initial directory:
    os.chdir(initdir)


def make_pdf(session, two_acq=False):
    """Method to make the PDF for the spider.

    :return: None
    """
    list_slices = [3, 6, 9, 12]
    slices = {'0': list_slices}
    labels = {'0': 'ADC Acquisition 1'}
    vmins = {'0': 0}
    vmaxs = {'0': 5 * 10**-9}
    nii = os.path.join(JOBDIR, '1', '%s_FIT_ADC_1.nii.gz' % session)
    images = [nii]
    if len(ACQ_FILES) > 1 and os.path.exists(ACQ_FILES[1]):
        images.append(os.path.join(
            JOBDIR, '2', '%s_FIT_ADC_2.nii.gz' % session))
        labels['1'] = 'ADC Acquisition 2'
        slices['1'] = list_slices
        vmins['1'] = 0
        vmaxs['1'] = 5 * 10**-9

    pdf_final = os.path.join(JOBDIR, '%s_ADC_report.pdf' % session)
    spiders.plot_images(pdf_final, 1, images, 'Compute ADC Verdict',
                        image_labels=labels, slices=slices,
                        vmins=vmins, vmaxs=vmaxs, orient='ax')


def write_dicom(pixel_array, filename, ds_ori,
                series_number, sop_id, nb_acq,
                project, subject, session):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_ori: original pydicom object of the pixel_array
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    :param nb_acq: acquisition number
    :param project: project name on XNAT
    :param subject: subject name on XNAT
    :param session: session name on XNAT
    """
    # Set the DICOM dataset
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = ds_ori.SOPInstanceUID
    file_meta.ImplementationClassUID = ds_ori.file_meta.ImplementationClassUID
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0" * 128)

    # Copy UID from dicom
    ds.InstanceCreatorUID = ds_ori.InstanceCreatorUID
    ds.SOPClassUID = ds_ori.SOPClassUID
    # ds.ReferencedStudySequence = ds_ori.ReferencedStudySequence
    ds.StudyInstanceUID = ds_ori.StudyInstanceUID
    ds.SeriesInstanceUID = ds_ori.SeriesInstanceUID

    # Other tags to set
    ds.SeriesNumber = series_number
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-', '')\
                                                   .replace(':', '')\
                                                   .replace('.', '')\
                                                   .replace(' ', '')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = 'ADC Map %d' % nb_acq
    ds.SeriesDescription = 'ADC Map %d' % nb_acq
    ds.PatientName = subject
    ds.PatientID = session
    # Set SeriesDate/ContentDate
    now = datetime.date.today()
    ds.SeriesDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.ContentDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.Modality = 'MR'
    ds.StudyDescription = project
    ds.AcquisitionNumber = 1
    ds.SamplesperPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.11'
    ds.NumberOfFrames = pixel_array.shape[2]
    ds.PixelRepresentation = 0
    ds.HighBit = 15
    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.SmallestImagePixelValue = pixel_array.min()
    ds.LargestImagePixelValue = pixel_array.max()
    ds.Columns = pixel_array.shape[0]
    ds.Rows = pixel_array.shape[1]

    # Organise the array:
    pixel_array2 = np.zeros((pixel_array.shape[0] * pixel_array.shape[2],
                             pixel_array.shape[1]))
    for i in range(pixel_array.shape[2]):
        pixel_array2[pixel_array.shape[0] * i:pixel_array.shape[0] * (i + 1),
                     :] = pixel_array[:, :, i]

    # Set the Image pixel array
    if pixel_array2.dtype != np.uint8:
        pixel_array2 = pixel_array2.astype(np.uint8)

    ds.PixelData = pixel_array2.tostring()
    # Save the image
    ds.save_as(filename)


def convert_nifti_2_dicoms(nifti_file, sour_obj, output_folder, nbacq,
                           project, subject, session):
    """Convert 3D niftis into DICOM files.

    :param nifti_file: path to the nifti file
    :param sour_obj: pydicom object for the source dicom
    :param output_folder: folder where the DICOM files will be saved
    :param nbacq: Acquisition number
    :param project: project name on XNAT
    :param subject: subject name on XNAT
    :param session: session name on XNAT
    :return: None
    """
    if not os.path.isfile(nifti_file):
        raise Exception("NIFTY file %s not found." % nifti_file)
    # Naming
    filename = 'ADC_MAP_%d_original' % nbacq

    # Edit Niftys
    f_img = nib.load(nifti_file)
    f_data = f_img.get_data()
    # Rotation 270
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    # Normalizing data:
    dmin = 0.0
    dmax = 5 * 10**-9
    f_data[f_data < dmin] = dmin
    f_data[f_data > dmax] = dmax

    # Scale numbers to uint8
    f_data = (f_data - dmin) * 255.0 / (dmax - dmin)

    # Subtract by setting all value to nan:
    # f_data[mask] = np.nan

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    series_number = 87011 + 11 * (nbacq - 1)
    sop_id = sour_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1]) + '.'

    # Write the dicom
    filename = os.path.join(output_folder, '%s.dcm' % filename)
    write_dicom(f_data, filename, sour_obj,
                series_number, sop_id, nbacq,
                project, subject, session)


if __name__ == '__main__':
    main()
