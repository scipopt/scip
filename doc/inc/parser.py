'''
Created on 06.01.2014

@author: bzfhende
'''
import re

sectiontag = "SECTION:"
questiontag = "QUESTION:"
answertag = "ANSWER:"
labeltag = "LABEL:"

def formatitem((question, answer, label)):
   return """
          array(
              'question'=>'%s',
              'answer'=>'%s',
              'label'=>'%s'
             )""" % (question, answer, label)

def formatsection((section, items)):
   print section
   return """
          array(
             'title'=>%s,
             'content'=>array(%s),
             )""" % (repr(section), ",\n".join(map(formatitem, items.__iter__())))

def formatallsections(sections, sectionitems):
   print sections
   return """
   $faq = array(
        %s
        );""" % (",\n".join([formatsection((section, sectionitems[section])) for section in sections]))

if __name__ == '__main__':
   sectionitems = {}
   sections = []
   items = []
   mode = None
   currquestion = curranswer = currlabel = None

   with open("faqtext.txt", 'r') as currentfile:
      for line in currentfile:
         notagline = line
         for tag in [sectiontag, questiontag, answertag, labeltag]:
            notagline = re.sub(r"%s" % (tag), "", notagline, flags=re.MULTILINE).strip()

         notagline = re.sub(r"'", r"\'", notagline)
         print notagline
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
            print "Keys:" , sectionitems.keys()
            f = open("faqdata.php", 'w')
            f.write("<?php \n%s\n ?>" % formatallsections(sections, sectionitems))
            f.close()
