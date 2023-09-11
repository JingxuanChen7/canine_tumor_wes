import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Map;

public class SummarizeBwaWgs {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String SamFile = args[0];
		String ID = args[1];
		String file_output = args[2];
		int read_length = Integer.parseInt(args[3]);
		String Cancer_Type = args[4];
		String Status = args[5];
		getSummaryBwa(SamFile, ID, file_output, read_length, Cancer_Type, Status);
	}

	public static void getSummaryBwa(String SamFile, String ID, String file_output, int read_length, String Cancer_Type,
			String Status) throws IOException {
		double total = 0;
		double unique = 0;// 3
		double duplicate = 0;// 3
		double Onemapped = 0;// 5,9
		double incorrect = 0;// 1
		double unmapped = 0;// 13

		try (FileReader fileReader = new FileReader(SamFile)) {
			try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
				String strCurrentLine_sam;
				while ((strCurrentLine_sam = bufferedReader.readLine()) != null) {

					if (strCurrentLine_sam.startsWith("@")) {
						;
					} else {

						String[] file_lst = strCurrentLine_sam.split("\t");
						total += 1;
						int status = Integer.parseInt(file_lst[1]) % 16;

						if (status == 5 || status == 9) {
							Onemapped += 1;
						} else if (status == 1) {
							incorrect += 1;
						} else if (status == 13) {
							unmapped += 1;
						} else if (status == 3) {
							if (file_lst[11].contains("XT:")) {
								String status2 = file_lst[11].split(":")[2];
								if (status2.equals("U") || status2.equals("M")) {
									unique += 1;
								} else if (status2.equals("R")) {
									duplicate += 1;
								}
							}
						}
					}
				}
				bufferedReader.close();
				
				long dog_total_length = (long) (Math.pow(10, 9) * 2.41);
				double total_mean_coverage = (read_length * total) / dog_total_length;
				double uniq_mapped_mean_coverage = (read_length * unique) / dog_total_length;
				double pairs = total / 2;

				String unique_rate = Double.toString(unique / total);
				String dup_rate = Double.toString(duplicate / total);
				String Onemap_rate = Double.toString(Onemapped / total);
				String incorrect_rate = Double.toString(incorrect / total);
				String unmapped_rate = Double.toString(unmapped / total);
				BufferedWriter bw = null;
				FileWriter fw = new FileWriter(file_output);
				bw = new BufferedWriter(fw);
				// bw.write("ID\tTotal_pairs\tUniquely_mapped_rate\tRepeatedly_mapped_rate\t1read_mapped_rate\tIncorrectly_mapped_rate\tUnmapped_rate\tUniquely_mapped\tRepeatedly_mapped\t1read_mapped\tIncorrectly_mapped\tUnmapped\tTotal_reads\n');
				bw.write("file_name\t" + "total_reads_pairs\t" + "Unique_mapped_rate\t" + "unmapped_rate\t"
						+ "total_uniq_number\t" + "total_mean_coverage\t" + "unique_mapped_mean_coverage\t"
						+ "cancer_type\t" + "status" + "\n");
				bw.write(ID + "\t" + pairs + "\t" + unique_rate + "\t" + unmapped_rate + "\t" + 
						Double.toString(unique) + "\t" + Double.toString(total_mean_coverage) + "\t"
						+ Double.toString(uniq_mapped_mean_coverage) + "\t" + Cancer_Type + "\t" + Status+"\n");

				bw.close();
			}

		}
	}

}   