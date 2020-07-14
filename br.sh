for i in $(seq 989)
do
echo $i
python build_network.py
python run_network.py
done
