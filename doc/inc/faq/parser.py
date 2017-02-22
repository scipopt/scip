helptext='''
Created on 06.01.2014

this parser generates php-readable data out of faqtext.txt.

The format is strict, but very easy: The FAQ is separated into SECTIONs,
each of which may contain an arbitrary number of items, which consist (in this order)
of a QUESTION, an alpha-numeric, unique LABEL of your choice and
an ANSWER:

Example:
SECTION: General Questions about SCIP
    QUESTION: What is SCIP?
    LABEL:whatisscip
    ANSWER:
    <p>
      SCIP is a constraint integer program solver. It can be used as a framework for <b>branch-cut-and-price</b>
      and contains all necessary plugins to serve as a <b>standalone solver</b> for MIP, MINLP, and PBO.
    </p>
    ...

The example has the right format to go:
   - We have four keywords for sections, questions, labels and answers, resp.
   - an item first mentions a QUESTION:, which then gets a LABEL:, and finally, an ANSWER:. Do not mix up this order.
   - Each keyword ends with a colon ':' (!)
   - The label is a one-liner
   - line break after the 'ANSWER:'-tag
   - Questions and answers should be simple HTML text elements, do not provide specific classes etc. therein.

To make faq-generation for both home use and web use possible, extra-keywords should be used to identify
relative locations to the documentation and file extensions:
   - use PATHTODOC to denote relative paths to the documentation.
   - use LINKEXT to denote file-extensions.

The parser automatically interpretes this.

@author: Gregor Hendel
@author Matthias Miltenberger
'''
import re
import sys

import argparse
parser = argparse.ArgumentParser(epilog=helptext)
parser.add_argument('--faqtextfile', type=str, default='faqtext.txt', help='faq text file')
parser.add_argument('--linkext', type=str, default='shtml', help='file extension for internal links')

sectiontag = "SECTION:"
questiontag = "QUESTION:"
answertag = "ANSWER:"
labeltag = "LABEL:"

def formatitem((question, answer, label)):
   '''
   returns a fully formatted php array containing all item information
   '''
   return """
          array(
              'question'=>'%s',
              'answer'=>'%s',
              'label'=>'%s'
             )""" % (question, answer, label)

def formatsection((section, items)):
   '''
   returns a fully formatted array to represent an entire section together with its items
   '''
   return """
          array(
             'title'=>%s,
             'content'=>array(%s),
             )""" % (repr(section), ",\n".join(map(formatitem, items.__iter__())))

def formatallsections(sections, sectionitems):
   '''
   formats a list of sections
   '''
   return """
   $faq = array(
        %s
        );""" % (",\n".join([formatsection((section, sectionitems[section])) for section in sections]))

if __name__ == '__main__':
   '''
   main execution method of script.

   reads an faqtext.txt and creates the corresponding faqdata.php php-readable data file
   '''
   sectionitems = {}
   sections = []
   items = []
   mode = None
   currquestion = curranswer = currlabel = None
   args = parser.parse_args()
   substitutions = {"LINKEXT":args.linkext,
                    "PATHTODOC": "."}
   faqtext = args.faqtextfile
   # open the faqtext file and iterate over all lines. change mode depending on last seen tag
   with open(faqtext, 'r') as currentfile:
      for line in currentfile:
         notagline = line
         for tag in [sectiontag, questiontag, answertag, labeltag]:
            notagline = re.sub(r"%s" % (tag), "", notagline, flags=re.MULTILINE).strip()

         notagline = re.sub(r"'", r"\'", notagline)

         # use the substitutions to get working links
         for key, val in substitutions.items():
             notagline = re.sub(key, val, notagline)
         # print notagline
         if sectiontag in line:

            if mode is not None:
               items.append((currquestion, curranswer, currlabel))
               sectionitems[sections[-1]] = items[:]

            sections.append(notagline)
            mode = None
            items = []

         elif questiontag in line:
            if mode == answertag:
               items.append((currquestion, curranswer, currlabel))
            currquestion = notagline
            mode = questiontag
         elif labeltag in line:
            currlabel = notagline
            mode = labeltag
         elif answertag in line:
            curranswer = None
            mode = answertag
         else:

            if mode == questiontag:
               currquestion += "\n%s" % notagline
            elif mode == answertag:
               if curranswer is None:
                  curranswer = "%s" % (notagline)
               else:
                  curranswer = "%s\n%s" % (curranswer, notagline)

      else:
         if mode is not None:
            items.append((currquestion, curranswer, currlabel))
            sectionitems[sections[-1]] = items[:]
            # print "Keys:" , sectionitems.keys()
            f = open("faqdata.php", 'w')
            f.write("<?php \n%s\n ?>" % formatallsections(sections, sectionitems))
            f.close()
