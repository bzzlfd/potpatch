
# depends on the result of `1.lda/1.bulk/2.scf`

ln -s ../2.scf/atom.config atom.config
ln -s ../2.scf/OUT.VR IN.VR

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



