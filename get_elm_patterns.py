# -*- coding: utf-8 -*-  
""" Grab ELM patterns based on
    http://elm.eu.org/browse.html
"""
import os

def get_elm_pattern(elm, link):
    """ Follow the link and grab the ELM regular expression """

    os.system('wget http://elm.eu.org/'
              + link)
    elm_file = elm + '.html'
    if os.path.exists(elm_file):
        with open(elm_file) as f:
            line = f.readline()
            while line.find('Pattern:') == -1:
                    line = f.readline()
            pattern = f.readline().split('>')[1].lstrip().split('<')[0].strip()
            
            print elm + '\t' + pattern.replace('…','...').upper()

def extract_pages():
    """ Return {} of ELM name to HTML link """

    elm2link = {}
    with open('browse.html') as f:
        for line in f:
            if line.find('elmPages/') != -1:
                if line.find('.html') != -1:
                    link = line.split("href=\"")[1].split("\"")[0]
                    elm = link.split('.')[0].split('/')[1]
                    elm2link[elm] = link
    return elm2link


os.system('wget http://elm.eu.org/browse.html')
elm2link = extract_pages()
for elm in elm2link:
    get_elm_pattern(elm, elm2link[elm])

# clean up
os.system('rm *html*')
