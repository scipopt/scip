#!/usr/bin/env python3

helptext='''
Created on 06.01.2014

this parser generates php-readable data and static html out of faqtext.txt.

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

The parser automatically interprets this.

The written faqdata.php can be transformed with php into static HTML by use of localfaq.php.
The written faq.inc produces the same output (without the need for php).

@author Gregor Hendel
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

faqinc_header='''<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
@media (prefers-color-scheme: light) {
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.answer {
    background-color: #fff;
}
}
.answer {
    padding-left:   1em;
}
</style>

'''

def formatitem(item):
   '''
   returns a fully formatted php array containing all item information
   '''
   (question, answer, label) = item
   return """
          array(
              'question'=>'%s',
              'answer'=>'%s',
              'label'=>'%s'
             )""" % (question, answer, label)

def formatsection(sec):
   '''
   returns a fully formatted array to represent an entire section together with its items
   '''
   (section, items) = sec
   sectionlabel = "faq_{}".format(re.sub(r'[^a-zA-Z]', '', section).lower())
   return """
          array(
             'title'=>%s,
             'label'=>%s,
             'content'=>array(%s),
             )""" % (repr(section),
                     repr(sectionlabel),
                     ",\n".join(map(formatitem, items.__iter__())))

def formatallsections(sections, sectionitems):
   '''
   formats a list of sections
   '''
   return """
   $faq = array(
        %s
        );""" % (",\n".join([formatsection((section, sectionitems[section])) for section in sections]))

def write_faqinc(f, sections, sectiontimes):
   '''
   write main part of faq.inc
   '''

   # table of contents
   for section in sections :
      sectionlabel = "faq_{}".format(re.sub(r'[^a-zA-Z]', '', section).lower())
      f.write('<h3><a class="reveal_faq" href="#%s"><span class="fa fa-caret-right"></span> %s</a></h3>' % (sectionlabel, section));
      f.write('<ol>')
      for item in sectionitems[section] :
         question = item[0]
         label = item[2]
         question = re.sub(r"\\'", r"'", question)
         f.write('  <li>\n')
         f.write('    <a class="reveal_faq" href="#%s">\n' % label)
         f.write('      %s' % question)
         f.write('    </a>\n')
         f.write('  </li>\n  ')
      f.write('</ol>\n')
   f.write('<hr />\n')

   for section in sections :
      sectionlabel = "faq_{}".format(re.sub(r'[^a-zA-Z]', '', section).lower())
      f.write('<h3 id="%s" class="anchor">' % sectionlabel)
      f.write('<span class="fa fa-caret-right"></span> %s<a href="#" class="pull-right"><span title="go to top" class="fa fa-caret-up"></span></a></h3><ol>' % section)
      for item in sectionitems[section] :
         question = item[0]
         answer = item[1]
         # this undoes the replacement done in notagline = re.sub(r"'", r"\'", notagline) below
         # if the PHP output is dropped, then we could clean this up
         question = re.sub(r"\\'", r"'", question)
         answer = re.sub(r"\\'", r"'", answer)
         label = item[2]
         f.write('  <li id="%s" class="anchor">\n' % label)
         f.write('    <h4>\n')
         f.write('      <a class="reveal_faq" href="#%s">%s</a>\n\n' % (label, question))
         f.write('      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>\n')
         f.write('    </h4>\n')
         f.write('    <div id="%s_ans">\n' % label)
         f.write('      %s' % answer)
         f.write('    </div>\n')
         f.write('  </li>\n  ')
      f.write('</ol>\n')


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
            f = open("faq.inc", 'w')
            f.write(faqinc_header)
            write_faqinc(f, sections, sectionitems)
            f.close()
