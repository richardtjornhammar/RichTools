/*
Copyright 2018 RICHARD TJÖRNHAMMAR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
// COMPILE COMMAND: g++ tree.cc -o tree
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>

#include <cassert>
#include <vector>
#include <queue>
#include <iomanip>

#include <stack>
//#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <time.h>

typedef union NodeInformation {
        double d;
        char c;
} NodeInformation;

class Node {
	public :
		Node() {} ;
		Node( std::string id ){ set_id(id); };
		Node( std::string id, std::vector<std::string> link_ids ) {
			id_=id; add_links(link_ids);
		}

		void set_id( std::string id ) { id_=id; };
		void add_label( std::string label ) { label_=label; };
		void add_description( std::string description ) { description_=description ; };
		void add_link( std::string l_id ) { links_.push_back(l_id); };
		void add_links( std::vector<std::string> edges )
		{
			links_.clear();
			for (std::vector<std::string>::iterator w = edges.begin() ; w != edges.end(); ++w)
			{
				links_.push_back(*w);
			}
		};

		std::string id( void ) {return id_;};
		std::string label( void ) {return label_;};
		std::vector<std::string> get_edges(void){ return(links_); };
		void show( void ) {
			std::cout << "NODE["<< id_ << ","<<label_<<"] - "<<description_<<"\nLINKS:\n";
			for ( auto l : links_ )
				std::cout << l << " ";
			std::cout << std::endl;
		};
		~Node(){};
	private :
		NodeInformation data_ ;
		std::string id_ ;
		std::string label_ ;
		std::string description_ ;
		int level_; // SHORTEST DISTANCE FROM ROOT NODE
		std::vector<double> metrics_; // VARIOUS NUM DATA CONTAINER BELONGING TO THE NODE
		std::vector<std::string> links_; // CONNECTION VIA HASH 
		std::vector< Node * >    edges_; // CONNECTION VIA POINTERS
};

class NodeGraph {
	public:
		NodeGraph() {};
		NodeGraph( Node n ) { add(n); };
		void add( Node n );
		void add( std::string id, std::string label ) {
			Node n(id); n.add_label(label);
			add(n);
		};
		void add( std::string id, std::string label, std::vector<std::string> link_ids ) 
		{
			Node n(id); n.add_label(label); n.add_links(link_ids);
			add(n);
		};
		Node get_node(std::string id) {
    			auto search = GraphMap_.find(id);
    			if(search == GraphMap_.end())
			{
				Node loose_node(""); id = "";
				loose_node.add_label(id);
				return( loose_node );
			}
			return( search->second );
		}
		void set_root_id(std::string rid){ root_id_=rid; };
		std::string get_root_id(void) { return root_id_; }
		Node get_root(){return(get_node(root_id_));};

		std::string DepthFirstSearch  ( void );
		std::string BreadthFirstSearch( void );

		void show(void);
		~NodeGraph() {};
	private:
		std::string root_id_;
		int Num_Edges_;
		int Num_Vertices_;
		std::unordered_map<std::string,Node> GraphMap_;
};

std::string
NodeGraph::DepthFirstSearch ( void )
{
	std::stack<std::string> S;
	std::string path("");
	std::unordered_set<std::string> visited;

	S.push( get_root_id() );
	while ( !S.empty() )
	{
		std::string v = S.top(); S.pop();
		std::unordered_set<std::string>::const_iterator got = visited.find(v);
		if ( got == visited.end() )
		{
			visited.insert(v);
			//path += v;
			Node ncurrent = get_node(v);
			std::vector<std::string> edges = ncurrent.get_edges();
			for (std::vector<std::string>::iterator w = edges.begin() ; w != edges.end(); ++w)
			{
				S.push(*w);
			}
			std::string l_label = ncurrent.label();
			path += l_label;
		}
	}
	return( path );
};

std::string
NodeGraph::BreadthFirstSearch ( void )
{
	std::queue<std::string> S;
	std::string path("");
	std::unordered_set<std::string> visited;

	S.push( get_root_id() );
	while ( !S.empty() )
	{
		std::string v = S.front(); S.pop();

		Node ncurrent = get_node(v);
		visited.insert(v);
		std::string l_label = ncurrent.label();
		path += l_label;
		std::vector<std::string> edges = ncurrent.get_edges();
		for (std::vector<std::string>::iterator w = edges.begin() ; w != edges.end(); ++w)
		{
			std::unordered_set<std::string>::const_iterator got = visited.find(*w);
			if ( got == visited.end() )
			{
				S.push(*w);
			}
		}

	}

	return( path );
};

void NodeGraph::add(Node n)
{
	assert(n.id().length()>0 && "Node is missing ID");
	GraphMap_.emplace (n.id(),n);
	if ( GraphMap_.size() == 1 )
	{
		set_root_id(n.id());
	}
}

void NodeGraph::show(void){

	for(auto kv : GraphMap_ )
	{
		std::cout << "\nMAPKEY:" << kv.first <<"\n";
		kv.second.show();
	}
}

/*
int main( int argc, char *argv[] )  // char
{
	int n_num=0 , n_chr=0 ;
	std::cout << "HELLO WORLD" << std::endl;
	for ( int i=1 ; i<argc ; i++ ) {
		std::cout << argv[i] << "\n";
	}
	std::string S,id,label;
	std::vector<std::string> v_ids = {"1","2"};
	S="R-HSA"; Node n( S ); S="HUMAN"; n.add_label( S );
	n.show();

	NodeGraph ng(n);

	S="R-HSA-0001"; Node n2( S ); S="HUMAN"; n2.add_label( S );
	S="R-HSA"; n2.add_link(S);
	ng.add(n2);
	ng.show( );

	NodeGraph RichTree;
	id="0"; label="2"; v_ids={"1","6"};
	RichTree.add(id,label,v_ids);
	id="1"; label="7"; v_ids={"2","3"};
	RichTree.add(id,label,v_ids);
	id="2"; label="2"; v_ids={"",""};
	RichTree.add(id,label,v_ids);
	id="3"; label="6"; v_ids={"4","5"};
	RichTree.add(id,label,v_ids);
	id="4"; label="5"; v_ids={"",""};
	RichTree.add(id,label,v_ids);
	id="5"; label="11"; v_ids={"",""};
	RichTree.add(id,label,v_ids);
	id="6"; label="5"; v_ids={"7",""};
	RichTree.add(id,label,v_ids);
	id="7"; label="9"; v_ids={"8",""};
	RichTree.add(id,label,v_ids);
	id="8"; label="4";
	RichTree.add(id,label,v_ids);
	id="9"; label="3"; v_ids={"",""};
	RichTree.add(id,label,v_ids);
	id="10"; label="1"; v_ids={"",""};
	RichTree.add(id,label,v_ids);

	id="6";
	Node nr = RichTree.get_node(id);
	nr.show();

	id="-6";
	Node nr2 = RichTree.get_node(id);
	nr2.show();

	std::cout << "ROOT::" << RichTree.get_root_id() << std::endl;

	std::string path = RichTree.DepthFirstSearch();
	std::cout << "ROUTE:: "<< path  << std::endl;

	std::string route = RichTree.BreadthFirstSearch();
	std::cout << "ROUTE:: "<< route << std::endl;

	return 0;
}
*/
