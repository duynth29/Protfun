import sys
import re
import os
import glob

numres = 1000  # number of residues per line
desclen = 1000  # maximum number of characters in nameline
ARGC = len(sys.argv)

if ARGC < 2:
    sys.exit("""
reformat.py from HHsuite3
Read a multiple alignment in one format and write it in another format
Usage: reformat.py [informat] [outformat] infile outfile [options] 
  or   reformat.py [informat] [outformat] 'fileglob' .ext [options] 

Available input formats:
   fas:     aligned fasta; lower and upper case equivalent, '.' and '-' equivalent
   a2m:     aligned fasta; inserts: lower case, matches: upper case, deletes: '-',
            gaps aligned to inserts: '.'
   a3m:     like a2m, but gaps aligned to inserts MAY be omitted
   sto:     Stockholm format; sequences in several blocks with sequence name at 
            beginning of line (hmmer output)
   psi:     format as read by PSI-BLAST using the -B option (like sto with -M first -r)
   clu:     Clustal format; sequences in several blocks with sequence name at beginning 
            of line
Available output formats:
   fas:     aligned fasta; all gaps '-'
   a2m:     aligned fasta; inserts: lower case, matches: upper case, deletes: '-', gaps 
            aligned to inserts: '.'
   a3m:     like a2m, but gaps aligned to inserts are omitted
   sto:     Stockholm format; sequences in just one block, one line per sequence
   psi:     format as read by PSI-BLAST using the -B option 
   clu:     Clustal format
If no input or output format is given the file extension is interpreted as format 
specification ('aln' as 'clu')

Options:
  -v int    verbose mode (0:off, 1:on)
  -num      add number prefix to sequence names: 'name', '1:name' '2:name' etc
  -noss     remove secondary structure sequences (beginning with >ss_)
  -sa       do not remove solvent accessibility sequences (beginning with >sa_)
  -M first  make all columns with residue in first sequence match columns 
            (default for output format a2m or a3m)
  -M int    make all columns with less than X% gaps match columns 
            (for output format a2m or a3m)
  -r        remove all lower case residues (insert states) 
            (AFTER -M option has been processed)
  -r int    remove all lower case columns with more than X% gaps
  -g ''     suppress all gaps
  -g '-'    write all gaps as '-'
  -uc       write all residues in upper case (AFTER all other options have been processed)
  -lc       write all residues in lower case (AFTER all other options have been processed)
  -l        number of residues per line (for Clustal, FASTA, A2M, A3M formats) 
            (default={})
  -d        maximum number of characers in nameline (default={})

Examples: reformat.py 1hjra.a3m 1hjra.a2m  
          (same as reformat.py a3m a2m 1hjra.a3m 1hjra.a2m)
          reformat.py test.a3m test.fas -num -r 90
          reformat.py fas sto '*.fasta' .stockholm
""".format(numres, desclen))

informat = ""
outformat = ""
infile = ""
outfile = ""
num = 0  # don't use sequence number as prefix: '>n|name'
noss = 0  # don't remove secondary structure
nosa = 1  # remove solvent accessibility sequences
options = " ".join(sys.argv[1:])
names = []  # names of sequences read in
seqs = []  # residues of sequences read in
n = 0  # number of sequences read in
k = 0  # counts sequences for output
remove_inserts = 0
remove_gapped = 0
matchmode = ""  # do not change capitalization
match_gaprule = 0  # columns with less than this percentage of gaps will be match columns
v = 2
update = 0
nss = -1  # index of secondary structure sequence
lname = None  # length of sequence name in clustal, stockholm, psi format
titleline = ""  # first line beginning with "#" in A3M, A2M, or FASTA files

informats = ["fas", "a2m", "a3m", "sto", "psi", "clu"]
outformats = ["fas", "a2m", "a3m", "sto", "psi", "clu", "ufas"]
gap = "default"
case = "default"

# Process options
if re.search(r" -i\s+(\S+) ", options):
    infile = re.search(r" -i\s+(\S+) ", options).group(1)
    options = re.sub(r" -i\s+\S+ ", " ", options)
if re.search(r" -o\s+(\S+) ", options):
    outfile = re.search(r" -o\s+(\S+) ", options).group(1)
    options = re.sub(r" -o\s+\S+ ", " ", options)
if re.search(r" -num ", options):
    num = 1
    desclen = 505
    options = re.sub(r" -num ", " ", options)
if re.search(r" -noss ", options):
    noss = 1
    options = re.sub(r" -noss ", " ", options)
if re.search(r" -sa ", options):
    nosa = 0
    options = re.sub(r" -sa ", " ", options)
if re.search(r" -g\s+'?(\S*)'? ", options):
    gap = re.search(r" -g\s+'?(\S*)'? ", options).group(1)
    options = re.sub(r" -g\s+'?\S*'? ", " ", options)
if re.search(r" -r\s+(\d+) ", options):
    remove_gapped = int(re.search(r" -r\s+(\d+) ", options).group(1))
    options = re.sub(r" -r\s+\d+ ", " ", options)
if re.search(r" -r ", options):
    remove_inserts = 1
    options = re.sub(r" -r ", " ", options)
if re.search(r" -lc ", options):
    case = "lc"
    options = re.sub(r" -lc ", " ", options)
if re.search(r" -uc ", options):
    case = "uc"
    options = re.sub(r" -uc ", " ", options)
if re.search(r" -v\s*(\d+) ", options):
    v = int(re.search(r" -v\s*(\d+) ", options).group(1))
    options = re.sub(r" -v\s*\d+ ", " ", options)
if re.search(r" -v ", options):
    v = 2
    options = re.sub(r" -v ", " ", options)
if re.search(r" -M\s+(\d+) ", options):
    matchmode = "gaprule"
    match_gaprule = int(re.search(r" -M\s+(\d+) ", options).group(1))
    options = re.sub(r" -M\s+\d+ ", " ", options)
if re.search(r" -M\s+first ", options):
    matchmode = "first"
    options = re.sub(r" -M\s+first ", " ", options)
if re.search(r" -u ", options):
    update = 1
    options = re.sub(r" -u ", " ", options)
if re.search(r" -l\s+(\S+) ", options):
    numres = int(re.search(r" -l\s+(\S+) ", options).group(1))
    options = re.sub(r" -l\s+\S+ ", " ", options)
if re.search(r" -lname\s+(\S+) ", options):
    lname = int(re.search(r" -lname\s+(\S+) ", options).group(1))
    options = re.sub(r" -lname\s+\S+ ", " ", options)
if re.search(r" -d\s+(\S+) ", options):
    desclen = int(re.search(r" -d\s+(\S+) ", options).group(1))
    options = re.sub(r" -d\s+\S+ ", " ", options)

# Assign informat, outformat, infile, and outfile
if not outfile:
    match = re.search(r"(\S+)\s*$", options)
    if match:
        outfile = match.group(1)
        options = re.sub(r"\S+\s*$", "", options)
    else:
        sys.exit("Error: no output file given: '{}'".format(options))

if not infile:
    match = re.search(r"(\S+)\s*$", options)
    if match:
        infile = match.group(1)
        options = re.sub(r"\S+\s*$", "", options)
    else:
        sys.exit("Error: no input file given: '{}'".format(options))

match = re.search(r"(\S+)\s*$", options)
if match:
    outformat = match.group(1)
    options = re.sub(r"\S+\s*$", "", options)
else:
    match = re.search(r"\S*\.(\S+?)$", outfile)
    if match:
        outformat = match.group(1).lower()
        if outformat in ["aln", "fa", "fasta", "afa", "afas", "afasta"]:
            outformat = "fas"
    else:
        print("Using FASTA output format: '{}'".format(options))
        outformat = "fas"

match = re.search(r"(\S+)\s*$", options)
if match:
    informat = match.group(1)
    options = re.sub(r"\S+\s*$", "", options)
else:
    match = re.search(r"\S*\.(\S+?)$", infile)
    if match:
        informat = match.group(1).lower()
        if informat in ["aln", "fa", "fasta"]:
            informat = "fas"
    else:
        print("Using FASTA input format: '{}'".format(options))
        informat = "fas"

# Warn if unknown options found
if options.strip():
    print("\nWARNING: unknown options '{}'".format(options.strip()))

# Check if input and output formats are valid
if informat not in informats:
    sys.exit("\nError: {} is not a valid input format option\n".format(informat))
if outformat not in outformats:
    sys.exit("\nError: {} is not a valid output format option\n".format(outformat))

if outformat == "ufas":
    gap = ""

def reformat(infile, outfile):
    global nss, titleline, n, names, seqs
    nss = -1
    titleline = ""

    # Input part
    with open(infile, "r") as f:
        lines = f.readlines()

    if informat in ["fas", "a2m", "a3m"]:
        n = 0
        seq = lines[0]
        if seq.startswith("#"):
            titleline = seq
            seq = seq.split("\n", 1)[1]
        seq = re.sub(r"\n\#.*", "", seq)

        if seq != ">":
            root_name = os.path.basename(infile).split(".")[0]
            names.append(root_name)
            seq = re.sub(r"[\n ]", "", seq)
            seqs.append(seq)
            n += 1

        for seq in lines[1:]:
            seq = re.sub(r"\n\#.*", "", seq)
            while ">" in seq:
                seq += next(lines)
            if re.match(r"^aa_", seq):
                continue
            if re.match(r"^sa_", seq) and nosa:
                continue
            if re.match(r"^ss_", seq):
                if noss:
                    continue
                nss = n
            seq = seq.strip()
            if seq:
                names.append(seq)
                seqs.append(seq)
                n += 1

    elif informat == "sto":
        seqhash = {}
        first_block = True
        n = 0

        for line in lines:
            line = line.strip()
            line = re.sub(r"\s+", " ", line)
            if re.match(r"^\#=GC SS_cons", line):
                line = re.sub(r"^\#=GC SS_cons", "ss_dssp", line)
            if re.match(r"^\#", line):
                continue
            if re.match(r"^\/\/", line):
                break
            if not line:
                first_block = False
                continue
            match = re.match(r"^\s*(\S+)\s+(\S+)", line)
            if not match:
                sys.exit("\nERROR found in stockholm format: {}".format(line))
            if match.group(1) not in seqhash:
                if re.match(r"^aa_", line):
                    continue
                if re.match(r"^sa_", line) and nosa:
                    continue
                if re.match(r"^ss_", line):
                    if noss:
                        continue
                    nss = n
                names.append(match.group(1))
                seqs.append(match.group(2))
                seqhash[match.group(1)] = n
                n += 1
                first_block = True
            else:
                if first_block:
                    sys.exit("\nERROR: sequence {} appears more than once per block\n".format(match.group(1)))
                seqs[seqhash[match.group(1)]] += match.group(2)

    elif informat == "clu":
        residues_per_line = 50
        block = 1
        n = 0
        k = 0

        for line in lines:
            line = line.strip()
            if re.match(r"CLUSTAL", line, re.IGNORECASE):
                continue
            if re.match(r"^\#", line):
                continue
            if re.match(r"^\/\/", line):
                break
            if not line:
                if k:
                    if n and n != k:
                        sys.exit("\nError: different number of sequences in blocks 1 and {} of {}\n".format(block, infile))
                    block += 1
                    n = k
                    k = 0
                continue
            match = re.match(r"^(\S+)\s+([ a-zA-Z0-9.-]+?)(\s+\d+)?$", line)
            if not match:
                if re.match(r"^[*.: ]*$", line):
                    continue
                if noss and (re.match(r"^aa_", line) or re.match(r"^ss_", line) or re.match(r"^sa_", line)):
                    continue
                line = line.strip()
                match = re.match(r"^(\S{1,20})([a-zA-Z0-9.-]{" + str(residues_per_line) + r"})(\s+\d+)?$", line)
                if not match:
                    sys.exit("\nError found in Clustal format in {}, line {}: '{}'\n".format(infile, line))
                name = match.group(1)
                residues = match.group(2)
                print("WARNING: Found no space between name and residues in {}, line {}: '{}'\n".format(infile, line))
            else:
                if noss and (re.match(r"^aa_", line) or re.match(r"^ss_", line) or re.match(r"^sa_", line)):
                    continue
                if re.match(r"^aa_", line) or re.match(r"^sa_", line):
                    continue
                if re.match(r"^ss_", line):
                    nss = n
                name = match.group(1)
                residues = match.group(2)
                residues = re.sub(r" ", "", residues)
                residues_per_line = len(residues)

            if block == 1:
                names.append(name)
                seqs.append(residues)
            else:
                seqs[k] += residues
                if names[k] != name:
                    print("WARNING: name of sequence {} in block 1 ({}) is not the same as in block {} ({}) in {}\n".format(k, names[k], block, name, infile))
            k += 1

        if k and n and n != k:
            sys.exit("\nError: different number of sequences in blocks 1 and {} of {}\n".format(block, infile))
        if not n:
            n = k

    elif informat == "psi":
        block = 1
        n = 0
        k = 0

        for line in lines:
            line = line.strip()
            if not line:
                if k:
                    if n and n != k:
                        sys.exit("\nError: different number of sequences in blocks 1 and {} of {}\n".format(block, infile))
                    block += 1
                    n = k
                    k = 0
                continue

            if noss and (re.match(r"^aa_", line) or re.match(r"^ss_", line) or re.match(r"^sa_", line)):
                continue
            if re.match(r"^aa_", line) or re.match(r"^sa_", line):
                continue
            if re.match(r"^ss_", line):
                nss = n
            match = re.match(r"^(\S+)\s+([ a-zA-Z0-9.-]+?)(\s+\d+)?$", line)
            name = match.group(1)
            residues = match.group(2)
            residues = re.sub(r" ", "", residues)

            if block == 1:
                names.append(name)
                seqs.append(residues)
            else:
                seqs[k] += residues
                if names[k] != name:
                    print("WARNING: name of sequence {} in block 1 ({}) is not the same as in block {} ({}) in {}\n".format(k, names[k], block, name, infile))
            k += 1

        if k and n and n != k:
            sys.exit("\nError: different number of sequences in blocks 1 and {} of {}\n".format(block, infile))
        if not n:
            n = k

    # Empty input file?
    if n == 0:
        sys.exit("\nERROR: input file {} contains no sequences\n".format(infile))

    # Transforming to upper-case
    if informat not in ["a3m", "a2m"]:
        seqs = [seq.upper() for seq in seqs]

    # Removing non-alphanumeric symbols
    seqs = [re.sub(r"[^A-Za-z0-9.~]", "", seq).replace("~", "-") for seq in seqs]

    # Filling up with gaps '.' or deleting gaps
    if informat == "a3m" and (not remove_inserts or matchmode):
        print("inserting gaps...")
        len_ins = []
        for seq in seqs:
            inserts = re.split(r"([A-Z]|-|~|[0-9])", "#" + seq + "#")
            j = 0
            for insert in inserts:
                if j >= len(len_ins):
                    len_ins.append(len(insert))
                elif len(insert) > len_ins[j]:
                    len_ins[j] = len(insert)
                j += 1

        for i, seq in enumerate(seqs):
            inserts = re.split(r"([A-Z]|-|~|[0-9])", "#" + seq + "#")
            j = 0
            for insert in inserts:
                inserts[j] = insert + "." * (len_ins[j] - len(insert))
                j += 1
            seqs[i] = "".join(inserts).replace("#", "")

    # Match state assignment
    if matchmode == "" and outformat in ["a3m", "a2m"]:
        matchmode = "first"

    if matchmode == "gaprule":
        gaps = []
        for seq in seqs:
            residues = list(seq)
            for l, residue in enumerate(residues):
                if residue in [".", "-"]:
                    if l >= len(gaps):
                        gaps.append(1)
                    else:
                        gaps[l] += 1

        for i, seq in enumerate(seqs):
            residues = list(seq)
            new_residues = ""
            for l, residue in enumerate(residues):
                if l >= len(gaps) or gaps[l] < 0.01 * match_gaprule * n:
                    if residue == ".":
                        new_residues += "-"
                    else:
                        new_residues += residue.upper()
                else:
                    if residue == "-":
                        new_residues += "."
                    else:
                        new_residues += residue.lower()
            seqs[i] = new_residues

    if matchmode == "first":
        match = []
        for k, name in enumerate(names):
            if not re.match(r"^(ss_|aa_|sa_)", name):
                break
        residues = list(seqs[k])
        for l, residue in enumerate(residues):
            if residue in [".", "-"]:
                match.append(0)
            else:
                match.append(1)

        for i, seq in enumerate(seqs):
            residues = list(seq)
            new_residues = ""
            for l, residue in enumerate(residues):
                if match[l]:
                    if residue == ".":
                        new_residues += "-"
                    else:
                        new_residues += residue.upper()
                else:
                    if residue == "-":
                        new_residues += "."
                    else:
                        new_residues += residue.lower()
            seqs[i] = new_residues

    # Remove gaps etc.
    if remove_gapped:
        gaps = []
        for seq in seqs:
            residues = list(seq)
            for l, residue in enumerate(residues):
                if residue in ["-", "."]:
                    if l >= len(gaps):
                        gaps.append(1)
                    else:
                        gaps[l] += 1

        for i, seq in enumerate(seqs):
            residues = list(seq)
            new_residues = ""
            for l, residue in enumerate(residues):
                if l >= len(gaps) or gaps[l] < 0.01 * remove_gapped * n:
                    new_residues += residue
            seqs[i] = new_residues

    seqs = [re.sub(r" ", "", seq) for seq in seqs]
    if remove_inserts:
        seqs = [re.sub(r"[a-z.]", "", seq) for seq in seqs]

    nin = n
    for i in range(n):
        if re.sub(r"[a-zA-Z0-9]", "", seqs[i]) == "":
            if v >= 2:
                print("Sequence contains only gaps and is removed: {}".format(names[i]))
            seqs.pop(i)
            names.pop(i)
            n -= 1

    names = [name[:desclen] for name in names]

    if outformat == "a3m":
        seqs = [re.sub(r"\.", "", seq) for seq in seqs]
    elif outformat in ["fas", "clu", "sto", "psi"]:
        seqs = [re.sub(r"\.", "-", seq) for seq in seqs]

    if gap != "default":
        seqs = [re.sub(r"\.", gap, seq).replace("-", gap) for seq in seqs]

    if case == "uc":
        seqs = [seq.upper() for seq in seqs]
    elif case == "lc":
        seqs = [seq.lower() for seq in seqs]

    # Check that sequences have same length
    if outformat not in ["a3m", "ufas"]:
        length = len(seqs[0])
        for i, seq in enumerate(seqs[1:], start=1):
            if len(seq) != length:
                print("\nError: Sequences in {} do not all have same length, e.g. >{}  (len={})  and  >{}  (len={})\n".format(
                    infile, names[0], length, names[i], len(seq)))
                if v >= 3:
                    print("{} {}\n{} {}\n\n".format(names[0], seqs[0], names[i], seq))
                sys.exit(1)

    # Remove html tags
    names = [re.sub(r"<[A-Za-z\/].*?>", "", name) for name in names]

    # Output part
    ndssp = -1
    nsa = -1
    npred = -1
    nconf = -1
    nquery = -1
    for i, name in enumerate(names):
        if re.match(r"^ss_dssp", name):
            ndssp = i
        elif re.match(r"^sa_dssp", name):
            nsa = i
        elif re.match(r"^ss_pred", name):
            npred = i
        elif re.match(r"^ss_conf", name):
            nconf = i
        elif nquery == -1 and not re.match(r"^aa_", name):
            nquery = i

    with open(outfile, "w") as f:
        if outformat in ["sto", "psi"]:
            if outformat == "sto":
                f.write("# STOCKHOLM 1.0\n\n")
                if not lname:
                    lname = 32
                if re.match(r"^\S+\s+(.*)", names[nquery]):
                    f.write("#=GF DE {}\n".format(re.match(r"^\S+\s+(.*)", names[nquery]).group(1)))
                refline = seqs[nquery]
                refline = re.sub(r"[a-z]", "-", refline)
                f.write("#=GC RF {}\n".format(refline))
                if ndssp >= 0:
                    f.write("#=GC SS_cons {}\n".format(seqs[ndssp]))

            if num:
                num = 2
                for i, name in enumerate(names):
                    if i in [ndssp, npred, nconf, nquery]:
                        continue
                    names[i] = re.sub(r"^(\S+)\#\d+", r"\1", name)
                    names[i] = re.sub(r"^(\S{1,25})\S+", r"\1#{}".format(num), name)
                    num += 1

            for i, name in enumerate(names):
                if i in [ndssp, npred, nconf]:
                    continue
                match = re.match(r"\s*(\S+)", name)
                if not lname:
                    lname = 32
                f.write("{:<{}} {}\n".format(match.group(1), lname, seqs[i]))

            if outformat == "sto":
                f.write("//\n")

        elif outformat == "clu":
            f.write("CLUSTAL\n\n\n")
            if num:
                num = 2
                for i, name in enumerate(names):
                    if i in [ndssp, npred, nconf, nquery]:
                        continue
                    names[i] = re.sub(r"^(\S+)\#\d+", r"\1", name)
                    names[i] = re.sub(r"^(\S{1,10})\S*", r"\1#{}".format(num), name)
                    num += 1

            while seqs[0]:
                for i, name in enumerate(names):
                    match = re.match(r"\s*(\S+)", name)
                    seqs[i] = re.sub(r"(\S{1," + str(numres) + r"})", r"\1\n", seqs[i])
                    if not lname:
                        lname = 18
                    f.write("{:<{}} {}\n".format(match.group(1), lname, seqs[i]))
                f.write("\n")

        else:
            if num:
                num = 2
                for i, name in enumerate(names):
                    if i in [ndssp, npred, nconf, nquery]:
                        continue
                    names[i] = re.sub(r"^(\S+)\#\d+", r"\1", name)
                    names[i] = re.sub(r"^(\S{1,25})\S+", r"\1#{}".format(num), name)
                    num += 1

            if titleline and outformat == "a3m":
                f.write("{}\n".format(titleline))

            for i, name in enumerate(names):
                seqs[i] = re.sub(r"(\S{" + str(numres) + r"})", r"\1\n", seqs[i])
                f.write(">{}\n{}\n".format(name, seqs[i]))

    if v >= 2:
        if nin == 1:
            print("Reformatted {} with 1 sequence from {} to {} and written to file {}".format(infile, informat, outformat, outfile))
        else:
            if nin != n:
                print("Removed {} sequences which contained no residues".format(nin - n))
            print("Reformatted {} with {} sequences from {} to {} and written to file {}".format(infile, n, informat, outformat, outfile))

if re.search(r"\*|\.\S*$", infile) or re.search(r"^\.\S*$", outfile):
    outext = re.search(r"\.(\S*)$", outfile).group(1)
    infiles = glob.glob(infile)
    print("{} files to reformat".format(len(infiles)))
    for infile in infiles:
        if not re.search(r"(\S+)\.\S+", infile):
            infile = re.search(r"(\S+)", infile).group(1)
        outfile = "{}.{}".format(re.search(r"(\S+)", infile).group(1), outext)
        if update and os.path.exists(outfile):
            continue
        if v >= 3:
            print("Reformatting {} from {} to {} ...".format(infile, informat, outformat))
        reformat(infile, outfile)
else:
    if v >= 3:
        print("Reformatting {} from {} to {} ...".format(infile, informat, outformat))
    reformat(infile, outfile)