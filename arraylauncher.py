from pathlib import Path
from time import sleep
import datetime
import re
import subprocess

USER = "cpetak"
SLEEP_TIME = 30
GUIDES_FOLDER = "guide_files"


def extract_chromosome(filename):
    pattern = r"guide_file_(NW_\d+\.\d+)\.txt"
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    return None


def get_time():
    now = datetime.datetime.now()
    return now.strftime("%Y-%m-%d %H:%M:%S")  # 2023-04-15 14:30:45

def run_check_output_safely(command):
    try:
        output = subprocess.check_output(command, stderr=subprocess.PIPE, text=True)
        return output
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(f"Error message: {e.stderr}")
        return None
    except Exception as e:
        print(f"Failed to execute command: {e}")
        return None

if __name__ == "__main__":

  guides = list(Path(GUIDES_FOLDER).glob("guide_file_*.txt"))
  print("Guides:", guides)
  print()

  while len(guides) > 0:
      guide = guides[0]
      guides = guides[1:]

      lines = guide.read_text().strip().split("\n")

      mychr = extract_chromosome(guide.as_posix())
      input_vcf = f"{mychr}_filtered.vcf_phased.vcf_processed.recode.vcf.gz"

      CMD = f"sbatch --array=1-{len(lines)} run_geva.sh guide_files/guide_file_{mychr}.txt 0.01 geva_results {input_vcf}"

      parts = CMD.split(" ")

      print(f"[{get_time()}] Processing {guide}")
      #print(input_vcf)
      print(f"[CMD] {CMD}")
      # print(parts)

      success = False
      wait = SLEEP_TIME
      while not success:
          output = run_check_output_safely(parts)
          if not output:
              print(
                  f"    [{get_time()}] Submission FAILED, sleeping for {wait} seconds"
              )
              sleep(wait)
              wait *= 2
          else:
              success = True
      print(f"    [{get_time()}] Submitted {guide}")
      print()

  print(f"[{get_time()}] DONE")
