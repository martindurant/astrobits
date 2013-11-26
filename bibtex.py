import urllib
import os

urlbase = """http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?db_key=AST&db_key=PRE&qform=AST&arxiv_sel=astro-ph&arxiv_sel=cond-mat&arxiv_sel=cs&arxiv_sel=gr-qc&arxiv_sel=hep-ex&arxiv_sel=hep-lat&arxiv_sel=hep-ph&arxiv_sel=hep-th&arxiv_sel=math&arxiv_sel=math-ph&arxiv_sel=nlin&arxiv_sel=nucl-ex&arxiv_sel=nucl-th&arxiv_sel=physics&arxiv_sel=quant-ph&arxiv_sel=q-bio&sim_query=YES&ned_query=YES&adsobj_query=YES&aut_logic=OR&obj_logic=OR&object=&start_mon=&end_mon=&ttl_logic=OR&title=&txt_logic=OR&text=&nr_to_return=200&start_nr=1&jou_pick=ALL&ref_stems=&data_and=ALL&group_and=ALL&start_entry_day=&start_entry_mon=&start_entry_year=&end_entry_day=&end_entry_mon=&end_entry_year=&min_score=&sort=SCORE&data_type=BIBTEX&aut_syn=YES&ttl_syn=YES&txt_syn=YES&aut_wt=1.0&obj_wt=1.0&ttl_wt=0.3&txt_wt=3.0&aut_wgt=YES&obj_wgt=YES&ttl_wgt=YES&txt_wgt=YES&ttl_sco=YES&txt_sco=YES&version=1"""

def parse_entry(textin):
    """Decompose a textual citation into a number of names and a year, for
turning into an ADS query"""
    text = textin.replace("\\","").replace(r"{","").replace(r"}","")
    bits = text.replace(","," ").split()
    name = ""
    names = "^"
    if textin=='' or textin[0]=='@':
        return
    for word in bits:
        if word.isalnum() and not(word.isalpha()):
            year = int(word)
            break
        if name == "": continue
        if name =="et": continue
        if name =="al": continue
        if name =="al.": continue
        if len(word)<=2:
            name = name+"%2C+"+word[0]
            names = names  + name+ "%0D%0A"
            name = ""
        elif name!="":
            names = names + name+ "%0D%0A" 
            name = word
        else:
            name = word
    return names, year

def get_ADS_bibtex(bibcode):
    """For the given bibcode key, get relavant bibtex entry"""
    url = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=%s&data_type=BIBTEX&db_key=AST"%bibcode
    print url
    temp_file,message = urllib.urlretrieve(url)
    fred = open(temp_file)
    for line in fred:
        print line,

def find_updated_key(bibcode):
    """If a reference has been updated in the database (e.g., arxiv->journal
publication), this enables you to find the updated key."""
    url = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=%s&data_type=BIBTEX&db_key=AST"%bibcode
    temp_file,message = urllib.urlretrieve(url)
    fred = open(temp_file)
    for line in fred:
        if line[:8] == "@ARTICLE":
            print 'Now known as ',line[9:-2]

import shutil
def remove_from_database(database,bibcode):
    """Backup bibtex (~) file, and make new version without bibcode included."""
    shutil.copy(database,database+'~')
    fred = open(database+'~')
    freda = open(database,'w')
    key=" "
    for line in fred:
        if "@" in line:
            key=line.split("{")[1][:-2]
        if key == bibcode:
            continue
        freda.write(line)
    fred.close()
    freda.close()

def add_to_database(database,bibcode,replace=False):
    """Get given bibcode from ADS, and add to the database, if it's not already present.
Replace involves first removing the key (and making a backup), and then
adding again - useful in case reference has changed. Note that in that
case the key might have changed too."""
    keys = get_keys(database)
    if bibcode in keys:
    	if not(replace):
	        return "Key already in database - nothing done"
	else: 
		print "Removing",bibcode
		remove_from_database(database,bibcode)
    url = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=%s&data_type=BIBTEX&db_key=AST"%bibcode
    temp_file,message = urllib.urlretrieve(url)
    fred = open(temp_file)
    freda = open(database,"a")
    for line in fred:
        if "@" in line:
            print "Adding",line.split("{")[1][:-2]
            break
    freda.write(line)
    for line in fred:
        freda.write(line)
    fred.close()
    freda.close()

def call_ADS(names,year):
    """Do an ADS search for the given author string (in HTML format) and year"""
    url = urlbase + "&author=" + names + "&start_year=%i&end_year=%i"%(year,year)
    temp_file,message = urllib.urlretrieve(url)
    return temp_file

def parse_listing(ads_output_file,output_stream):
    """Take a file given by ADS, and trim it for the first bibtex entry"""
    fred = open(ads_output_file)
    for line in fred:
        if "@" in line:
            output_stream.write(line)
            break
    for line in fred:
        output_stream.write(line)
        if line == "}\n":
            output_stream.write("\n")
            break
    fred.close()

def alt_find_bibitems(tex_file,marker = r"\bibitem"):
    r"""As find_bibitem, but for freetext form. Each reference has "marker" in it."""
    fred = open(tex_file)
    bibitems=[]
    start = ""
    for line in fred:
        if marker in line:
            break
    for line in fred:
        if line=='\n' or line[0]=='%': continue
        if not(marker in line):
            start = start + line[:-1]
        else:
            bibitems.append(start)
            start=""
    bibitems.append(start)
    return bibitems

def find_bibitems(tex_file,marker = r"\bibitem"):
    r"""Scan through a LaTeX file, extracting the main string
from bibliographic definitions. This for \bibitem[arg]{short}{full}
structure"""
    fred = open(tex_file)
    bibitems=[]
    for line in fred:
        if not(marker in line):
            continue
        start = line.split("{")[-1]
        if not("}\n" in line): 
            for line in fred:
                start = start.strip("\n") + line.strip("\n")
                if "}" in line:
                    break
        bibitems.append(start)
    fred.close()
    return bibitems

def find_bibitems3(tex_file,marker = r"\bibitem"):
    r"""Scan through a LaTeX file, extracting the main string
from bibliographic definitions. This for \bibitem[arg]{short}{full}
structure"""
    fred = open(tex_file)
    bibitems=[]
    for line in fred:
        if not(marker in line):
            continue
        start = line.split("}")[-1].strip("\n")
        bibitems.append(start)
    fred.close()
    return bibitems

def get_keys(bibfile):
    """find all the ADS bib keys in the given file"""
    fred = open(bibfile)
    keys = []
    for line in fred:
        if "@" in line:
            keys.append(line.split("{")[1][:-2])
    return keys

def merge_from_tex(texfile,bibfile,style=0,marker = r'\bibitem'):
    """Take the natbib definitions in texfile and merge them into the bibtex
file given. Style choses between two ways to do natbibing"""
    fred = open(bibfile,"a")
    keys = get_keys(bibfile)
    if style==0:
        bibitems = find_bibitems(texfile,marker)
    elif style==1:
        bibitems = alt_find_bibitems(texfile,marker)
    else:
        bibitems = find_bibitems3(texfile,marker)
    for item in bibitems:
        try:
            out = parse_entry(item)
            if out:
                names,year = out
            else: continue
            print names,year
            temp_file = call_ADS(names,year)
            key = get_keys(temp_file)
            if key[0] in keys:
                print "Skipping %s"%key[0]
                continue
            parse_listing(temp_file,fred)
            os.remove(temp_file)
        except KeyboardInterrupt:
            break
        except:
            print "Error in ", item
            pass
    fred.close()



