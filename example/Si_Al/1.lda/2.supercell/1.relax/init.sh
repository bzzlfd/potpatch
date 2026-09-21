
potpatch mksupcl -i ../../1.bulk/2.scf/atom.config -o atom.config -r 1.0 -s 4 4 4
mv atom.config _atom_.config
awk '
$1==14 && $2==0 && $3==0 && $4==0 && !done {
  sub(/14/, "13")
  done=1
}
{ print }
' _atom_.config > atom.config
rm _atom_.config
