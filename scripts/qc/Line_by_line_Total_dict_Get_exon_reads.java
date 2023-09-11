import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Collections;

// The java version of the target_reads.py scripts.
// The script include BWA analysis result
// Each sample takes around 15 min vs 2 hr in python script
public class Line_by_line_Total_dict_Get_exon_reads {

	public static void main(String[] args) throws IOException {

		String ExonInterval = "/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
		String SamFile = args[0];
		String file_name = args[1];
		File file_output = new File(args[2]);
		int read_length = Integer.parseInt(args[3]);
		String Cancer_Type = args[4];
		String Status = args[5];

		// Create Total CDS-Interval-Dictionary
		Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict = createCDSIntervalDict(ExonInterval);
		Set<String> Total_chrom_list = Total_interval_dict.keySet();

		for (String ele : Total_chrom_list) {
			Collections.sort(Total_interval_dict.get(ele));

		}

		ArrayList<String> summary = getExonCount(SamFile, read_length, Total_interval_dict);

		writeOutput(file_name, file_output, Cancer_Type, Status, summary);
	}

	/**
	 * @param file_name
	 * @param file_output
	 * @param Cancer_Type
	 * @param Status
	 * @param summary
	 * @throws IOException
	 */
	private static void writeOutput(String file_name, File file_output, String Cancer_Type, String Status,
			ArrayList<String> summary) throws IOException {
		BufferedWriter bw = null;
		FileWriter fw = new FileWriter(file_output);
		bw = new BufferedWriter(fw);

		bw.write("File_name\t" + "Total_Pairs\t" + "Uniq_Means\t" + "Uniq_mapped_rate\t"
				+ "Uniq_Exonic_region_mapped_rate\t" + "CancerType\t" + "Status\t" + "UnmappedRate\t"
				+ "DuplicateMapped_rate\t" + "Onemapped_rate\t" + "IncorrectMapped_rate\t" + "Total_line\t"
				+ "Total_unique\t" + "Total_pass\t" + "Total_Unmapped\t" + "Total_Duplicate\t" + "Total_Onemapped\t"
				+ "Total_Incorrect\t" + '\n');

		bw.write(file_name + '\t' + summary.get(0) + '\t' + summary.get(1) + '\t' + summary.get(2) + '\t'
				+ summary.get(3) + '\t' + Cancer_Type + '\t' + Status
				+ '\t'+ summary.get(4) + '\t'  + summary.get(5) + '\t' + summary.get(6) + '\t' + summary.get(7) + '\t' + summary.get(8) + '\t' + summary.get(9) + '\t'
				+ summary.get(10) + '\t' + summary.get(11) + '\t' + summary.get(12) + '\t' + summary.get(13) + '\t'
				+ summary.get(14) + '\n');
		bw.close();
	}

	/**
	 * @param ExonInterval
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws NumberFormatException
	 */
	private static Map<String, ArrayList<ExonLocationInfo>> createCDSIntervalDict(String ExonInterval)
			throws FileNotFoundException, IOException, NumberFormatException {
		Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict = new HashMap<>();
		BufferedReader objReader = null;
		objReader = new BufferedReader(new FileReader(ExonInterval));
		String strCurrentLine;
		while ((strCurrentLine = objReader.readLine()) != null) {
			ArrayList<ExonLocationInfo> ExonLocation = new ArrayList<ExonLocationInfo>();

			String[] tokens = strCurrentLine.split(":");
			String chrom = tokens[0];
			int start = Integer.parseInt(tokens[1].split("-")[0]);
			int end = Integer.parseInt(tokens[1].split("-")[1]);

			ExonLocation.add(new ExonLocationInfo(start, end));
			if (Total_interval_dict.containsKey(chrom) == false) {

				Total_interval_dict.put(chrom, ExonLocation);
			} else {
				Total_interval_dict.get(chrom).add(new ExonLocationInfo(start, end));

			}

		}
		objReader.close();
		return Total_interval_dict;
	}

	public static ArrayList<String> getExonCount(String SamFile, int read_length,
			Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict) throws NumberFormatException, IOException {
		int total_unique = 0;// 3 is the remainder of 16 in the flag
		int total_cds_pass = 0;
		int total_line = 0;
		int total_Duplicate = 0;// flag%16 = 3
		int total_Onemapped = 0;// flag%16 = 5,9
		int total_Incorrect = 0;// flag%16 =1
		int total_Unmapped = 0;// flag%16 =13
		int total_cds_pos = 35683639;
		long dog_total_length = (long) (Math.pow(10, 9) * 2.41);

		ArrayList<String> summary_list = new ArrayList<String>();
		try (FileReader fileReader = new FileReader(SamFile)) {
			try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {

				String strCurrentLine_sam;

				while ((strCurrentLine_sam = bufferedReader.readLine()) != null) {
					int[] summary = analyzeLine(strCurrentLine_sam, read_length, Total_interval_dict);
					// bufferedReader.close();
					total_line += summary[0];
					total_unique += summary[1];
					total_cds_pass += summary[2];
					total_Unmapped += summary[3];
					total_Duplicate += summary[4];
					total_Onemapped += summary[5];
					total_Incorrect += summary[6];

				}
				bufferedReader.close();
			}
		}
		
		double uniq_mapped_rate = Double.valueOf(total_unique) / Double.valueOf(total_line);
		double uniq_Exonic_region_mapped_rate = Double.valueOf(total_cds_pass) / Double.valueOf(total_unique);
		double unmappedRate = Double.valueOf(total_Unmapped) / Double.valueOf(total_line);
		double DuplicateMapped_rate = Double.valueOf(total_Duplicate) / Double.valueOf(total_line);
		double Onemapped_rate = Double.valueOf(total_Onemapped) / Double.valueOf(total_line);
		double IncorrectMapped_rate = Double.valueOf(total_Incorrect) / Double.valueOf(total_line);
		double uniq_mean_coverage = Double.valueOf(read_length) * Double.valueOf(total_cds_pass) / Double.valueOf(dog_total_length);

		int pairs = total_line / 2;
		summary_list.add(String.valueOf(pairs));
		summary_list.add(String.valueOf(uniq_mean_coverage));
		summary_list.add(String.valueOf(uniq_mapped_rate));
		summary_list.add(String.valueOf(uniq_Exonic_region_mapped_rate));
		summary_list.add(String.valueOf(unmappedRate));
		summary_list.add(String.valueOf(DuplicateMapped_rate));
		summary_list.add(String.valueOf(Onemapped_rate));
		summary_list.add(String.valueOf(IncorrectMapped_rate));
		summary_list.add(String.valueOf(total_line));
		summary_list.add(String.valueOf(total_unique));
		summary_list.add(String.valueOf(total_cds_pass));
		summary_list.add(String.valueOf(total_Unmapped));
		summary_list.add(String.valueOf(total_Duplicate));
		summary_list.add(String.valueOf(total_Onemapped));
		summary_list.add(String.valueOf(total_Incorrect));

		return summary_list;

	}

	// Analyzed for each line
	public static int[] analyzeLine(String strCurrentLine_sam, int read_length,
			Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict) {

		int numberOfline = 0;
		int cds_pass_line = 0;
		int unique = 0;
		int duplicate = 0;// flag%16 = 3
		int onemapped = 0;// flag%16 = 5,9
		int incorrect = 0;// flag%16 =1
		int unmapped = 0;// flag%16 =13
		int[] summary = new int[7];
		ArrayList<String> transcript_list = new ArrayList<String>();
		if (strCurrentLine_sam.startsWith("@")) {
			;
		} else {
			String[] file_lst = strCurrentLine_sam.split("\t");
			numberOfline += 1;
			// String reads_name = file_lst[0];
			String reads_chr = file_lst[2];
			int reads_position = Integer.parseInt(file_lst[3]);
			int status = Integer.parseInt(file_lst[1]) % 16;
			if (status == 3) {

				if (file_lst[11].contains("XT:")) {
					String status2 = file_lst[11].split(":")[2];
					if (status2.equals("U") || status2.equals("M")) {
						unique += 1;

						ArrayList<ExonLocationInfo> Total_exome_loc = Total_interval_dict.get(reads_chr);

						int location_status = binarySearch(Total_exome_loc, 0, (Total_exome_loc.size() - 1),
								reads_position, read_length);

						if (location_status != -2) {
							// transcript_list.add(reads_name);
							cds_pass_line += 1;
						}

					}

					else if (status2.equals("R")) {
						duplicate += 1;
					}
				}
			}

			else if (status == 5 || status == 9) {
				onemapped += 1;
			} else if (status == 1) {
				incorrect += 1;
			}

			else if (status == 13) {
				unmapped += 1;
			}

			summary[0] = numberOfline;
			summary[1] = unique;
			summary[2] = cds_pass_line;
			summary[3] = unmapped;
			summary[4] = duplicate;
			summary[5] = onemapped;
			summary[6] = incorrect;

		}

		return summary;

	}

	public static int withinRegion(int reads_position, int exom_start, int exom_end, int read_length) {
		int reads_end = reads_position + read_length - 1;
		int reads_start = reads_position;
		int value = 0;

		// means the read is greater
		if (reads_start > exom_end) {
			value = 0;
		}
		// means the read is smaller
		else if (reads_end < exom_start) {
			value = -1;
		}
		// means the read is in the region
		else {
			value = 1;
		}
		return value;

	}

	public static int binarySearch(ArrayList<ExonLocationInfo> arr, int left, int right, int reads_position,
			int read_length) {

		// Check base case
		if (right >= left) {

			int mid = (left + right) / 2;
			int exon_start_value = arr.get(mid).getstart();
			int exon_end_value = arr.get(mid).getend();
			// If element is present at the middle itself
			int region_status = withinRegion(reads_position, exon_start_value, exon_end_value, read_length);
			if (region_status == 1) {
				return mid;
			}
			// If element is smaller than mid, then it can only be present in left subarray
			else if (region_status == -1) {
				return binarySearch(arr, left, mid - 1, reads_position, read_length);
			}

			// Else the element can only be present in right subarray
			else if (region_status == 0) {
				return binarySearch(arr, mid + 1, right, reads_position, read_length);
			}
		}

		// Element is not present in the array
		return -2;

	}
}

class ExonLocationInfo implements Comparable<ExonLocationInfo> {

	int start;
	int end;

	public ExonLocationInfo(int start, int end) {
		this.start = start;
		this.end = end;
	}

	public int getstart() {
		return start;
	}

	public int getend() {
		return end;
	}

	public int compareTo(ExonLocationInfo exon) {
		// use start point to sort
		return this.start - exon.start;
	}

}
