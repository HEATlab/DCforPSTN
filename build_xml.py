import glob
import os


## \file build_xml.py
#  \brief convert an AMPL format model file to xml file


##
# \fn build_xml(xml_name: str, model_name: str):
# \brief convert an AMPL format model file to xml file
#
# @param xml_name       The filename of the output xml file
# @param model_name     The filename of the input model file
#
# @post An xml file generated from the model file that can be submitted to
#       the neos server
def build_xml(xml_name: str, model_name: str):
   # Read in the optimiziation problem specification
   model_file = open(model_name, 'r')
   model = model_file.read()
   model_file.close()

   # Build the XML document
   doc_beg = "<document>"
   doc_end = "</document>"
   choice = ("<category>minco</category>\n"
           "<solver>BARON</solver>\n"
           "<inputMethod>AMPL</inputMethod>")
   empty_data = "<data><![CDATA[]]></data>"
   empty_commands = "<commands><![CDATA[]]></commands>"
   empty_comments = "<comments><![CDATA[]]></comments>"

   template = (f"{doc_beg}\n{choice}\n\n"
           "<model><![CDATA[\n"
           f"{model}\n"
           "]]></model>\n\n"
           f"{empty_data}\n\n"
           f"{empty_commands}\n\n"
           f"{empty_comments}\n\n"
           f"{doc_end}")

   xml_file = open(xml_name, 'w')
   xml_file.write(template)
   xml_file.close()


##
# A short interactive program for converting files.
def main():
   path = input("Please enter the path of directory:\n")
   out_path = input("Please enter the path of out directory:\n")
   listOfFile = glob.glob(os.path.join(path, '*.mod'))

   for fname in listOfFile:
       p, f = os.path.split(fname)
       outname = os.path.join(out_path, f[:-4] + '.xml')
       print("Processing: ", f)
       build_xml(outname, fname)

if __name__ == "__main__":
    main()
