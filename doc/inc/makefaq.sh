#! /bin/sh

echo "creating faq.inc..."
python parser.py && php index.php > faq.inc
rm faqdata.php
echo "done"
