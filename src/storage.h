#ifndef _STORAGE_H_
#define _STORAGE_H_


#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;



class storage
{

public:

  void write_to_file(bool cycles);

  void do_averaging();
  
  void add_to_switches(unordered_map<string, int> switches);

  void add_to_time(map<int, int>& phen, map<int,int>& adult);

  map<string, int>& get_switch_tally();

  map<int, int>& get_phenotype_tally();

  map<int, int>& get_adult_tally();

  void add_to_edges(map<pair<int,int>,int>& tally);

  map<pair<int,int>,int>& get_edges();


  ~storage();


private:  

  map<string, int> switch_tally{};


  map<int, int> phenotype_tally{};

  map<int,int> adult_tally{};



  map<pair<int,int>, int> edge_tally{};




};


#endif