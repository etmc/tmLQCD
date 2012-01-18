#! /bin/bash
# SED script to change from home-rolled _Complex double to libc _Complex double
for var in "$@"
do
  cp "$var" "$var".bak

  # Taking care of the type 'complex'
    echo "Replace the header inclusion"
    egrep -n --color '\#include\s*\"complex.h\"' < "$var"
    sed -i"" -r 's/\#include\s*\"complex.h\"/\#include <complex.h>/g' "$var"
    
    echo "Replace all occurences of the type 'complex"
    egrep -n --color '([^<]|^)\<complex\>' < "$var"
    sed -i"" -r 's/([^<]|^)\<complex\>/\1_Complex double/g' "$var"

  # Taking care of all the '.re' and '.im' mess
    echo "Replace \"(*x).\" by \"x->\""
    egrep -n --color '\(\*(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\)\.' < "$var"
    sed -i"" -r 's/\(\*(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\)\./\1->/g' "$var"

    echo "Replace all complex equivalencies"
    # NB We need to use multiline searching here, so this one's a b*tch -- the arcance stuff is building this multiline buffer
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([\*\+\-\\]?=) *(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *\5\6\.im\;' < "$var"
    sed -i"" -r -n '1h;1!H;${;g;s/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([\*\+\-\\]?=) *(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *\5\6\.im\;/\1 \4 \5\6;/g;p;}' "$var"

    echo "Replace all complex equivalencies with conjugation"
    # NB We need to use multiline searching here, so this one's a b*tch -- the arcance stuff is building this multiline buffer
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *\5\6\.im\;' < "$var"
    sed -i"" -r -n '1h;1!H;${;g;s/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *\5\6\.im\;/\1 \4 \5 (\6);/g;p;}' "$var"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *-\5\6\.im\;' < "$var"
    sed -i"" -r -n '1h;1!H;${;g;s/([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*([A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *-(.*)(\<[A-Za-z_][A-Za-z0-9_]*\>(\[[^]]+\])?)\.re\;[ 	\r\n]*\1\.im *\4 *-\5\6\.im\;/\1 \4 -\5conj(\6);/g;p;}' "$var"

    echo "Replace all '.re' expressions with self-reference"
    egrep -n --color '((([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^\]]+\])?)\.re) *([\*\+\-\\]?=)[^\;]*\1[^\;]*\;' < "$var"
    sed -i"" -r 's/((([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^\]]+\])?)\.re) *([\*\+\-\\]?=)([^\;])*\1([^\;]*)\;/\2 \5 (\6creal(\2)\7) + cimag(\2) * I\;/g' "$var"

    echo "Replace all '.im' expressions with self-reference"
    egrep -n --color '((([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^\]]+\])?)\.im) *([\*\+\-\\]?=)([^\;])*\1([^\;])*\;' < "$var"
    sed -i"" -r 's/((([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^\]]+\])?)\.im) *([\*\+\-\\]?=)([^\;])*\1([^\;]*)\;/\2 \5 creal(\2) + (\6cimag(\2)\7) * I\;/g' "$var"

    echo "Replace all complex expressions with double assignment"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *([^\;]+)\;[ 	\r\n]*\1\.im *\4 *([^\;]+)\;' < "$var"
    sed -i"" -r -n '1h;1!H;${;g;s/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *([^\;]+)\;[	 \r\n]*\1\.im *\4 *([^\;]+)\;/\1 \4 (\5) + (\6) * I;/g;p;}' "$var"

    echo "Replace all '.re' expressions with assignment"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *([^\;]+)\;' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re *([-\+\*\/]?=) *([^\;]+)\;/\1 \4 (\5) + cimag(\1) * I;/g' "$var"
   
    echo "Replace all '.im' expressions with assignment"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.im *([-\+\*\/]?=) *([^\;]+)\;' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.im *([-\+\*\/]?=) *([^\;]+)\;/\1 \4 creal(\1) + (\5) * I;/g' "$var"

    echo "Replace all '.re' expressions without assignments"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re( *[^=]|==)' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.re( *[^=]|==)/creal(\1)\4/g' "$var"

    # Specific case for call of component on a function / macro output
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*\([^\)]*\))\.re( *[^=]|==)' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*\([^\)]*\))\.re( *[^=]|==)/creal(\1)\3/g' "$var"

    echo "Replace all '.im' expressions without assignments"
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.im( *[^=]|==)' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)\.im( *[^=]|==)/cimag(\1)\4/g' "$var"

    # Specific case for call of component on a function / macro output
    egrep -n --color '(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*\([^\)]*\))\.im( *[^=]|==)' < "$var"
    sed -i"" -r 's/(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*\([^\)]*\))\.im( *[^=]|==)/creal(\1)\3/g' "$var"
  # Replacing the macros
    echo "Replace calls to _add_assign_complex_conj(x,z,w)"
    egrep -n --color '_add_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_add_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) += conj(\2) * (\3)/g' "$var"

    echo "Replace calls to _add_assign_complex(x,z,w)"
    egrep -n --color '_add_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_add_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) += (\2) * (\3)/g' "$var"

    echo "Replace calls to _diff_assign_complex_conj(x,z,w)"
    egrep -n --color '_diff_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_diff_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) -= conj(\2) * (\3)/g' "$var"

    echo "Replace calls to _diff_assign_complex(x,z,w)"
    egrep -n --color '_diff_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_diff_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) -= (\2) * (\3)/g' "$var"

    echo "Replace calls to _mult_assign_complex_conj(x,z,w)"
    egrep -n --color '_mult_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_mult_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = conj(\2) * (\3)/g' "$var"

    echo "Replace calls to _minus_mult_assign_complex_conj(x,z,w)"
    egrep -n --color '_minus_mult_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_minus_mult_assign_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = -conj(\2) * (\3)/g' "$var"

    echo "Replace calls to _mult_assign_complex(x,z,w)"
    egrep -n --color '_mult_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_mult_assign_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2) * (\3)/g' "$var"

    echo "Replace calls to _add_mult_complex(x,z,w)"
    egrep -n --color '_add_mult_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_add_mult_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) += (\2) * (\3)/g' "$var"

    echo "Replace calls to _diff_complex(z,w)"
    egrep -n --color '_diff_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_diff_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) -= (\2)/g' "$var"

    echo "Replace calls to _add_complex(z,w)"
    egrep -n --color '_add_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_add_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) += (\2)/g' "$var"

    echo "Replace calls to _mult_real(z,w,r)"
    egrep -n --color '_mult_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_mult_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2) * (\3)/g' "$var"

    echo "Replace calls to _mult_assign_real(z,r)"
    egrep -n --color '_mult_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_mult_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) *= (\2)/g' "$var"

    echo "Replace calls to _complex_square_norm(z)"
    egrep -n --color '_complex_square_norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_square_norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)/creal(\1 * conj(\1))/g' "$var"

    echo "Replace calls to _complex_norm(z)"
    egrep -n --color '_complex_norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_square_norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)/cabs(\1)/g' "$var"

    echo "Replace calls to _complex_zero(z)"
    egrep -n --color '_complex_zero[ 	]*\([ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_zero[ 	]*\([ 	]*([^,\)]+)[ 	]*\)/(\1) = 0.0/g' "$var"

    echo "Replace calls to _complex_one(z)"
    egrep -n --color '_complex_one[ 	]*\([ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_one[ 	]*\([ 	]*([^,\)]+)[ 	]*\)/(\1) = 1.0/g' "$var"

    echo "Replace calls to _complex_1norm(z)"
    egrep -n --color '_complex_1norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_square_1norm[ 	]*\([ 	]*([^,\)]+)[ 	]*\)/(fabs(creal(\1)) + fabs(cimag(\1)))/g' "$var"

    echo "Replace calls to _complex_set(z,r,s)"
    egrep -n --color '_complex_set[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_set[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2) + (\3) * I/g' "$var"

    echo "Replace calls to _div_complex(x,z,w)"
    egrep -n --color '_div_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's#_div_complex[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)#(\1) = (\2) / (\3)#g' "$var"

    echo "Replace calls to _div_real(z,w,r)"
    egrep -n --color '_div_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's#_div_real[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)#(\1) = (\2) / (\3)#g' "$var"

    echo "Replace calls to _complex_conj(z,w)"
    egrep -n --color '_complex_conj[ 	]*\(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])?)([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = conj(\2)/g' "$var"

    echo "Replace calls to _complex_chgsig(z,w)"
    egrep -n --color '_complex_chgsig[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_chgsig[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = -(\2)/g' "$var"

    echo "Replace calls to _assign(z,w)"
    egrep -n --color '_assign[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_assign[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2)/g' "$var"

    echo "Replace calls to _complex_conj_chgsig(z,w)"
    egrep -n --color '_complex_conj_chgsig[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_complex_conj_chgsig[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = -conj(\2)/g' "$var"

    echo "Replace calls to _assign_add_conj(z,v,w)"
    egrep -n --color '_assign_add_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_assign_add_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2) + conj(\3)/g' "$var"

    echo "Replace calls to _assign_diff_conj(z,v,w)"
    egrep -n --color '_assign_diff_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)' < "$var"
    sed -i"" -r 's/_assign_diff_conj[ 	]*\([ 	]*([^,\)]+),[ 	]*([^,\)]+),[ 	]*([^,\)]+)[ 	]*\)/(\1) = (\2) - conj(\3)/g' "$var"

#    echo "Clean up parentheses around atomic expressions"
#    egrep -n --color '(^|[ 	]+)(-?)\([ 	]*(-?[ABD-Za-z_][0-9A-Za-z_]*|C[0-9A-Za-z_]+|-?[0-9]+\.?[0-9]*)[ 	]*\)' < "$var"
#    sed -i"" -r 's/(^|[ 	]+)(-?)\([ 	]*(-?[ABD-Za-z_][0-9A-Za-z_]*|C[0-9A-Za-z_]|-?[0-9]+\.?[0-9]**)[ 	]*\)/\1\2\3/g' "$var"

    echo "Remove redundant imaginary components"
    egrep -n --color '[ 	]*[\+\-][ 	]*\(?0(\.0*)?\)?[ 	]*\*[ 	]*I' < "$var"
    sed -i"" -r 's/[ 	]*[\+\-][ 	]*\(?0(\.0*)?\)?[ 	]*\*[ 	]*I//g' "$var"

    echo "Remove multiplication by inverse"
    egrep -n --color '\*[ 	]*\(1\./(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])*(\.re|\.im)?)\)' < "$var"
    sed -i"" -r 's#\*[ 	]*\(1\./(([A-Za-z_][A-Za-z0-9_]*\.|[A-Za-z_][A-Za-z0-9_]*\->)*[A-Za-z_][A-Za-z0-9_]*(\[[^]]+\])*(\.re|\.im)?)\)#/ \1#g' "$var"

    echo "Clean up parentheses in assignment"
    egrep -n --color '([^=])=[ 	]*\(([^\(\)]+)\)[ 	]*\;' < "$var"
    sed -i"" -r 's/([^=])=[ 	]*\(([^\(\)]+)\)[ 	]*\;/\1= \2;/g' "$var"

    echo "Clean up awkward spacing in parentheses"
    egrep -n --color '([\(\)])[ 	]+([\(\)])' < "$var"
    sed -i"" -r 's/([\(\)])[ 	]+([\(\)])/\1\2/g' "$var"
done

