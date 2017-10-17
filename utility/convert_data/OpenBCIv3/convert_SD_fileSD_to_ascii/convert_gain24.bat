for %%f in (*.TXT) do gawk -f "convert_gain24.awk" "%%f" > "%%f.new.csv"


