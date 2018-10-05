layout_network: main_layout_network.cpp ForceLayout.cpp Quadtree.cpp
	g++ -std=c++11 -o layout_network $^

clean:
	rm layout_network

clean_layouts:
	rm *.layout

%.layout: ./layout_network %.edges
	$^ > $@

plot4x4.png: drawlayout.R ~/Downloads/sample_networks
	Rscript $^ $@

obsplot-%.png: obsmodel.R ~/Downloads/sample_networks
	Rscript $^ $@

.PHONY: obsplots

obsplots: obsplot-85.png obsplot-19156.png

# TODO: more flexibility NxM, paging via pdf call, make vectorized?
# TODO: Roz will want filtering condition - accomplish via separate script