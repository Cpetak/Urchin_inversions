from bs4 import BeautifulSoup
import sys
import csv


def extract_intergenic_value(filename):
    # Read the HTML file
    with open(filename, "r") as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, "html.parser")

    # Find all td elements with intergenic_region text
    for td in soup.find_all("td"):
        if td.b and td.b.text.strip() == "intergenic_region":
            # Get the next td with class='numeric'
            numeric_td = td.find_next("td", class_="numeric")
            if numeric_td:
                # Remove comma and return as integer
                return int(numeric_td.text.strip().replace(",", ""))

    return None


if len(sys.argv) < 2:
    print("Please provide the input HTML file as an argument")
    sys.exit(1)

filename = sys.argv[1]
outfilename = sys.argv[2]
intergenic_value = extract_intergenic_value(filename)
#print(intergenic_value)

with open(outfilename, mode='a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow([intergenic_value])

