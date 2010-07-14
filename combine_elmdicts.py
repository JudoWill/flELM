"""When simplifying the residue alphabet,
   I made so many ELM patterns that I had
   to split up the file. This combines the
   results.

   Enter results directory & maximum suffix """
import sys, os, global_settings

dir = sys.argv[1] # working/Jul12/
end = int(sys.argv[2]) # 10

for g in global_settings.TEST_GENOMES:
    cat_line = ''
    for x in xrange(end+1):
        file = dir + 'elmdict_' + g + '.init' + str(x) + ' '
        cat_line += file
    os.system('cat ' + cat_line + '> '
              + dir + 'elmdict_' + g + '.init')
