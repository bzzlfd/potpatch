
# depends on the result of `2.hse/1.bulk/2.scf`

ln -s ../2.scf/atom.config atom.config
ln -s ../2.scf/OUT.VR IN.VR
for i in `ls ../2.scf/OUT.HSE*`; do
  file=`basename $i`
  ln -s ../2.scf/$file $file
done

cat > gen.kpt << EOF
commet line
10
0.0 0.0 0.0 Gamma
0.5 0.0 0.0 X
10
0.5 0.0 0.0 X
0.5 0.5 0.0 M
EOF
split_kp.x gen.kpt



