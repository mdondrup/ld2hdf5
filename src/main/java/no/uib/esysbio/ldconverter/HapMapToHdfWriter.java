package no.uib.esysbio.ldconverter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Collection;
import java.util.Iterator;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;
import org.apache.commons.io.*;

/**
 * This program converts a hapmap file <a
 * href="ftp://ftp.ncbi.nlm.nih.gov/hapmap/ld_data/latest/00README.txt">Go here
 * for more explanation</a>
 * 
 */
public class HapMapToHdfWriter {

	static char delimiter = ' ';

	/** Name for the dataset containing the snp id */

	private final static String DNAME_SNP = "/snp.id.one";

	/** Name for the dataset containing the corresponding snp. */
	private final static String DNAME_COR_SNP = "/snp.id.two";

	private final static String DNAME_VALUES = "/value.set";
	private final static String DNAME_MARKER_ONE = "/marker.pos.one";
	private final static String DNAME_MARKER_TWO = "/marker.pos.two";
	private final static int RANK = 1;
	private final static long[] MAX_DIMS = { HDF5Constants.H5S_UNLIMITED };

	/*
	 * Creates a new dataset to add values in.
	 */
	private static void createIntegerDataset(int fid, long[] dims,
			long[] chunkSize, String groupName) throws Exception {
		int did_snp = -1, did_CorSNP = -1, did_value = -1, did_marker_one_pos = -1, did_marker_two_pos = -1, type_int_id = -1, type_double_id = -1, sid = -1, plist = -1, group_id = -1;

		try {

			type_int_id = H5.H5Tcopy(HDF5Constants.H5T_STD_I32LE);
			type_double_id = H5.H5Tcopy(HDF5Constants.H5T_IEEE_F64LE);
			// use variable length to save space

			sid = H5.H5Screate_simple(RANK, dims, MAX_DIMS);

			// figure out creation properties

			plist = H5.H5Pcreate(HDF5Constants.H5P_DATASET_CREATE);

			H5.H5Pset_layout(plist, HDF5Constants.H5D_CHUNKED);

			H5.H5Pset_chunk(plist, RANK, chunkSize);

			H5.H5Pset_deflate(plist, 6);

			group_id = H5.H5Gcreate(fid, groupName, HDF5Constants.H5P_DEFAULT);

			did_snp = H5.H5Dcreate(group_id, groupName + DNAME_SNP,
					type_int_id, sid, plist);

			did_marker_one_pos = H5.H5Dcreate(group_id, groupName
					+ DNAME_MARKER_ONE, type_int_id, sid, plist);

			did_marker_two_pos = H5.H5Dcreate(group_id, groupName
					+ DNAME_MARKER_TWO, type_int_id, sid, plist);

			did_CorSNP = H5.H5Dcreate(group_id, groupName + DNAME_COR_SNP,
					type_int_id, sid, plist);

			did_value = H5.H5Dcreate(group_id, groupName + DNAME_VALUES,
					type_double_id, sid, plist);

			System.out.println("created for chr " + groupName);
			// H5.H5Pset_szip(plist, HDF5Constants.H5_SZIP_NN_OPTION_MASK, 8);
		} finally {
			try {
				H5.H5Pclose(plist);

			} catch (HDF5Exception ex) {
			}
			try {
				H5.H5Sclose(sid);

			} catch (HDF5Exception ex) {
			}
			try {
				H5.H5Dclose(did_snp);
				H5.H5Dclose(did_CorSNP);
				H5.H5Dclose(did_value);
				H5.H5Dclose(did_marker_one_pos);
				H5.H5Dclose(did_marker_two_pos);
				H5.H5Gclose(group_id);
			} catch (HDF5Exception ex) {
			}

		}
	}

	/**
	 * Operations for reading from a directory of text file, with suffix txt,
	 * parsing them and creating a corresponding HDF file for each of the files.
	 * 
	 * @param fid
	 *            File id created by hdf.
	 * @param sourceFolder
	 *            Folder to look for the text files in.
	 * @throws Exception
	 */
	private static void writeDataFromFileToInt(int fid, String sourceFolder)
			throws Exception {
		int did_SNP = -1, did_CorSNP = -1, did_value = -1, did_marker_one = -1, did_marker_two = -1, type_int_id, type_double_id = -1, msid = -1, fsid = -1, timesWritten = 0, group_id = -1;

		Collection<File> fileCollection = FileUtils.listFiles(new File(
				sourceFolder), new String[] { "txt" }, false);

		int filesAdded = 0;

		for (File sourceFile : fileCollection) {
			try {
				String chromosome_tmp = sourceFile.getName().split("_")[1];
				String chromosome = chromosome_tmp.substring(3,
						chromosome_tmp.length());
				chromosome = "/" + chromosome + "/";

				BigFile tmpfile = new BigFile(sourceFile.getAbsolutePath());

				int numLines = 0;
				long t1 = System.currentTimeMillis();
				for (String ldLine : tmpfile)
					numLines++;

				long t2 = System.currentTimeMillis();

				System.out.println("lines " + numLines + " time spent "
						+ (t2 - t1));
				long[] DIMS = { numLines };
				long[] CHUNK_SIZE = { 1000000 };
				int BLOCK_SIZE = 3000000;

				long[] count = { BLOCK_SIZE };
				createIntegerDataset(fid, DIMS, CHUNK_SIZE, chromosome);

				group_id = H5.H5Gopen(fid, chromosome);

				did_SNP = H5.H5Dopen(group_id, chromosome + DNAME_SNP);
				did_CorSNP = H5.H5Dopen(group_id, chromosome + DNAME_COR_SNP);
				did_value = H5.H5Dopen(group_id, chromosome + DNAME_VALUES);
				did_marker_one = H5.H5Dopen(group_id, chromosome
						+ DNAME_MARKER_ONE);
				did_marker_two = H5.H5Dopen(group_id, chromosome
						+ DNAME_MARKER_TWO);
				type_int_id = H5.H5Dget_type(did_SNP);
				type_double_id = H5.H5Dget_type(did_value);
				fsid = H5.H5Dget_space(did_SNP);
				msid = H5.H5Screate_simple(RANK, count, null);

				int[] currentSNPIdArray = new int[BLOCK_SIZE];
				int[] correspondingSNPIdArray = new int[BLOCK_SIZE];
				int[] markerOnePositionArray = new int[BLOCK_SIZE];
				int[] markerTwoPositionArray = new int[BLOCK_SIZE];
				double[] valueArray = new double[BLOCK_SIZE];

				BigFile ldFile = new BigFile(sourceFile.getAbsolutePath());

				int idx = 0, block_indx = 0, start_idx = 0;
				long t0 = 0;
				t0 = System.currentTimeMillis();
				System.out.println("Started to parse the file");

				int it = (int) (DIMS[0] / BLOCK_SIZE);
				int rest = (int) (DIMS[0] % BLOCK_SIZE);
				boolean lastWrite = false;
				int currentLine = 0;
				timesWritten = 0;
				for (String ldLine : ldFile) {
					char[] ldArray = new char[ldLine.length() - 1];
					ldLine.getChars(0, ldArray.length - 1, ldArray, 0);
					int numWhiteSpaces = 0;
					int stopOne = 0;
					int stopTwo = 0;

					int snpIdStart = 0;
					int snpIdEnd = 0;

					int snpIdTwoStart = 0;
					int snpIdTwoEnd = 0;

					int posOneStop = 0;

					int posTwoStop = 0;

					for (int j = 0; j < ldArray.length; j++) {
						if (ldArray[j] == delimiter)
							numWhiteSpaces++;

						if (numWhiteSpaces == 6 && stopOne == 0)
							stopOne = j;

						if (numWhiteSpaces == 7 && stopTwo == 0)
							stopTwo = j;

						if (numWhiteSpaces == 3 && snpIdStart == 0)
							snpIdStart = j;

						if (numWhiteSpaces == 4 && snpIdEnd == 0) {
							snpIdEnd = j;
						}

						if (numWhiteSpaces == 4 && snpIdTwoStart == 0)
							snpIdTwoStart = j;

						if (numWhiteSpaces == 5 && snpIdTwoEnd == 0)
							snpIdTwoEnd = j;

						if (numWhiteSpaces == 1 && posOneStop == 0)
							posOneStop = j;

						if (numWhiteSpaces == 2 && posTwoStop == 0)
							posTwoStop = j;
					}

					currentSNPIdArray[idx] = Integer.valueOf(ldLine.substring(
							snpIdStart + 3, snpIdEnd));

					correspondingSNPIdArray[idx] = Integer.valueOf(ldLine
							.substring(snpIdTwoStart + 3, snpIdTwoEnd));
					valueArray[idx] = Double.valueOf(ldLine.substring(
							stopOne + 1, stopTwo));

					markerOnePositionArray[idx] = Integer.valueOf(ldLine
							.substring(0, posOneStop));

					markerTwoPositionArray[idx] = Integer.valueOf(ldLine
							.substring(posOneStop + 1, posTwoStop));

					idx++;

					// Writes the data to the hdf file.
					if (idx == BLOCK_SIZE) {
						idx = 0;
						H5.H5Sselect_hyperslab(fsid,
								HDF5Constants.H5S_SELECT_SET,
								new long[] { start_idx }, null, count, null);
						H5.H5Dwrite(did_SNP, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, currentSNPIdArray);

						H5.H5Dwrite(did_CorSNP, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT,
								correspondingSNPIdArray);

						H5.H5Dwrite(did_marker_one, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT,
								markerOnePositionArray);

						H5.H5Dwrite(did_marker_two, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT,
								markerTwoPositionArray);

						H5.H5Dwrite(did_value, type_double_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, valueArray);

						block_indx++;
						start_idx = currentLine + 1;
						timesWritten++;
						if (timesWritten == it)
							lastWrite = true;
					}

					if (lastWrite && rest != 0 && idx == rest) {
						System.out.println("entries left " + idx
								+ " start idx " + start_idx);
						int[] tmp_SNP = new int[rest];
						int[] tmp_CORR = new int[rest];
						int[] tmp_marker_one = new int[rest];
						int[] tmp_marker_two = new int[rest];
						double[] tmp_values = new double[rest];
						for (int j = 0; j < rest; j++) {
							tmp_SNP[j] = currentSNPIdArray[j];
							tmp_CORR[j] = correspondingSNPIdArray[j];
							tmp_values[j] = valueArray[j];
							tmp_marker_one[j] = markerOnePositionArray[j];
							tmp_marker_two[j] = markerTwoPositionArray[j];
						}

						long[] lastDim = { rest };
						msid = H5.H5Screate_simple(RANK, lastDim, null);

						H5.H5Sselect_hyperslab(fsid,
								HDF5Constants.H5S_SELECT_SET,
								new long[] { start_idx }, null, lastDim, null);

						H5.H5Dwrite(did_SNP, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, tmp_SNP);

						H5.H5Dwrite(did_CorSNP, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, tmp_CORR);

						H5.H5Dwrite(did_marker_one, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, tmp_marker_one);

						H5.H5Dwrite(did_marker_two, type_int_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, tmp_marker_two);

						H5.H5Dwrite(did_value, type_double_id, msid, fsid,
								HDF5Constants.H5P_DEFAULT, tmp_values);
						System.out.println("start " + start_idx);
						block_indx++;

					}
					currentLine++;

				}
				filesAdded++;

				System.out.println("Finished parsing the file "
						+ (System.currentTimeMillis() - t0) + " for number "
						+ filesAdded);
				// FileUtils.deleteQuietly(sourceFile);
			} finally {
				try {
					H5.H5Gclose(group_id);
					H5.H5Sclose(fsid);

				} catch (HDF5Exception ex) {
				}
				try {
					H5.H5Sclose(msid);
				} catch (HDF5Exception ex) {
				}
				try {
					H5.H5Dclose(did_SNP);
					H5.H5Dclose(did_CorSNP);
					H5.H5Dclose(did_value);
					H5.H5Dclose(did_marker_one);
					H5.H5Dclose(did_marker_two);
				} catch (HDF5Exception ex) {
				}
			}
		}

	}

	public void parseHapMapFileToHdf(String sourceFolder, String resultFile)
			throws Exception {
		int fid = -1;

		// Creates hdf5 file, gives back a file handler used to refer to the
		// file later
		fid = H5.H5Fcreate(resultFile, HDF5Constants.H5F_ACC_TRUNC,
				HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

		if (fid < 0)
			return;

		try {
			File sourceFolderObject = new File(sourceFolder);

			if (!sourceFolderObject.exists())
				System.out.print("Folder " + sourceFolder + " does not exist");
			else
				writeDataFromFileToInt(fid, sourceFolder);
		} finally {
			H5.H5Fclose(fid);
		}
	}

	public static void main(String[] args) {
		try {
			/* Path to folder for hap map files. */
			String pathToFolder = args[0];
			/* Name of the file to create. */
			String hdfFileToCreate = args[1];
			HapMapToHdfWriter instance = new HapMapToHdfWriter();
			instance.parseHapMapFileToHdf(pathToFolder, hdfFileToCreate);
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

}
