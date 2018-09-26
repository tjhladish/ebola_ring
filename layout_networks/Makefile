layout_network:
	g++ -std=c++11 -o layout_network main_layout_network.cpp ForceLayout.cpp Quadtree.cpp

clean:
	rm layout_network

clean_layouts:
	rm *.layout

%.layout: ./layout_network %.edges
	$^ > $@

plot4x4.png: drawlayout.R ~/Downloads/sample_networks
	Rscript $^ $@