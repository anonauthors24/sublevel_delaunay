/* Copyright 2021 TU Graz
   Author: Michael Kerber
   
   This file is part of mpp_utils
   
   mpp_utils is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   mpp_utils is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public License
   along with mpp_utils.  If not, see <https://www.gnu.org/licenses/>.*/

#pragma once

#include<phat/boundary_matrix.h>

#include<algorithm>
#include<unordered_map>
#include <string>

#include <mpp_utils/basic.h>
#include <mpp_utils/Coordinate_traits_with_map.h>

#include<mpp_utils/Pre_column_struct.h>

namespace mpp_utils {

  template<typename Coordinate_>
  struct Grade_struct {
    typedef Coordinate_ Coordinate;
    Coordinate first_val;
    index first_index;
    Coordinate second_val;
    index second_index;
    Grade_struct() {}
    Grade_struct(Coordinate x, Coordinate y) : first_val(x), second_val(y) {}
    Grade_struct(index ind_x, index ind_y, Coordinate val_x, Coordinate val_y) : first_val(val_x), first_index(ind_x), second_val(val_y), second_index(ind_y) {}
    Grade_struct(const Grade_struct<Coordinate>& other) : first_val(other.first_val), first_index(other.first_index),
							  second_val(other.second_val), second_index(other.second_index) {}
    bool operator== (const Grade_struct<Coordinate>& other) {
      return this->first_val==other.first_val && this->second_val==other.second_val;
    }
  };

  template<typename Representation_=phat::vector_vector, 
           typename CoordinateTraits = Coordinate_traits_with_map<double> >
    class Graded_matrix  : public phat::boundary_matrix<Representation_> {
    
    public:

    typedef Representation_ Representation;
    typedef CoordinateTraits Coordinate_traits;
    
    typedef typename Coordinate_traits::Coordinate Coordinate;

    typedef Grade_struct<Coordinate> Grade;

    std::vector<Coordinate> x_vals,y_vals;

    index num_grades_x;
    index num_grades_y;

    index num_rows;
    
    std::vector<Grade> grades;

    std::vector<Grade> row_grades;

    bool grade_indices_assigned;

    /*
    Graded_matrix() {}

    Graded_matrix(Graded_matrix& other) {
      std::copy(other.x_vals.begin(),other.x_vals.end(),std::back_inserter(this->x_vals));
      std::copy(other.y_vals.begin(),other.y_vals.end(),std::back_inserter(this->y_vals));
      // not complete
    }
    */
    
    void print(bool print_row_grades=false,bool print_indices=true) {
      std::cout << "Number of columns: " << this->get_num_cols() << std::endl;
      std::cout << "Number of rows: " << this->num_rows << std::endl;

      for(int i=0;i<this->get_num_cols();i++) {
	if(print_indices) {
	  std::cout << grades[i].first_index <<  " " << grades[i].second_index << " -- ";
	} else {
	  std::cout << grades[i].first_val <<  " " << grades[i].second_val << " -- ";
	}
	phat::column col;
	this->get_col(i,col);
	for(std::size_t j=0;j<col.size();j++) {
	  std::cout << col[j] << " ";
	}
	std::cout << std::endl;

      }

      if(print_row_grades) {
	std::cout << "Row grades:" << std::endl;
	for(index j=0;j<this->num_rows;j++) {
	  std::cout << j << " : ";
	  if(print_indices) {
	    std::cout << row_grades[j].first_index <<  " " << row_grades[j].second_index << std::endl;
	  } else {
	    std::cout << row_grades[j].first_val <<  " " << row_grades[j].second_val << std::endl;
	  }
	}
      }
      
    }

    bool is_local(index i) {
      index p = this->get_max_index(i);
      //std::cout << "Info: " << i << " "<<  p  << "rows " << this->row_grades.size() << std::endl;
      return p!=-1 && (this->grades[i].first_val == this->row_grades[p].first_val) &&  (this->grades[i].second_val == this->row_grades[p].second_val) ;
      
    }

    bool pivot_is_dominating(index i) {
      if(this->is_empty(i)) {
	return false;
      }
      if(this->is_local(i)) {
	return true;
      }
      index p = this->get_max_index(i);
      Grade& pgr = this->row_grades[p];
      //std::cout << "#### NEW COLUMN ####" << std::endl;
      //std::cout << "pivot grade is " << pgr.first_index << ", " << pgr.second_index << std::endl;
      std::vector<index> col;
      this->get_col(i,col);
      for(int j=0;j<col.size();j++) {
	Grade& curr=this->row_grades[col[j]];
	//std::cout << "next grade is " << curr.first_index << ", " << curr.second_index << std::endl;
	if(curr.first_val>pgr.first_val || curr.second_val>pgr.second_val) {
	  //std::cout << "NOT HERE" << std::endl;
	  return false;
	}
      }
      
      return true;
    }


    template<typename OutStream>
    void print_in_rivet_format(OutStream& out,bool header=true,bool print_rows=true) {
      
      if(header) {
	out << "firep\nfirst parameter\nsecond parameter\n";

	out << this->get_num_cols() << " " << this->num_rows << " 0" << std::endl;
      }	

      Coordinate_traits& traits = Coordinate_traits::get_instance();

      for(index i=0;i<this->get_num_cols();i++) {
	//std::cout << "Printing " << grades[i].first_val << " - Result: " << traits.to_string(grades[i].first_val) << std::endl;
	//std::cout << grades[i].first_val <<  " " << grades[i].second_val << std::endl;
	traits.to_stream(out,grades[i].first_val);
	out <<  " ";
	traits.to_stream(out,grades[i].second_val);
	out <<  " ; ";
	std::vector<index> col;
	this->get_col(i,col);
	for(index j=0;j<col.size();j++) {
	  out << col[j] << " ";
	}
	out << "\n";
      }
      if(print_rows) {
	//std::cout << "Printing rows " << num_rows << " " << row_grades.size() << std::endl;
	for(index i=0;i<num_rows;i++) {
	  traits.to_stream(out,row_grades[i].first_val);
	  out << " ";
	  traits.to_stream(out,row_grades[i].second_val);
	  out << " ; \n";
	}
      }
      out << std::flush;
    }

  }; // of class Graded_matrix


  template<typename GradedMatrix>
    void check_boundaries(GradedMatrix& M, std::string msg) {
    for(index i=0;i<M.get_num_cols();i++) {
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=1;j<col.size();j++) {
	if(col[j-1]>=col[j]) {
	  std::cout << "Bad boundary " << col[j-1] << " " << col[j] << " at column " << i << "(" << msg << " matrix)" << std::endl;
	  std::exit(1);
	}
      }
    }
  }

  template<typename GradedMatrix>
    void check_grade_sanity(GradedMatrix& M) {

    typedef typename GradedMatrix::Grade Grade;
    
    for(index i=0;i<M.get_num_cols();i++) {
      Grade& col_gr = M.grades[i];
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	Grade& row_gr = M.row_grades[col[j]];
	if(row_gr.first_val>col_gr.first_val || row_gr.second_val>col_gr.second_val) {
	  std::cout << "Bad grading at: " << col[j] << "( " << row_gr.first_val << ", " << row_gr.second_val << ") " << i << "( " << col_gr.first_val << ", " << col_gr.second_val << ")" << std::endl;
	  std::exit(1);
	}
	assert(row_gr.first_val<=col_gr.first_val);
	assert(row_gr.second_val<=col_gr.second_val);
      }
    }
  }

  template<typename InputIterator>
    void copy_grades_to_rows(InputIterator begin,
			     InputIterator end) {

    typedef typename InputIterator::value_type GradedMatrix;
    auto it1=begin, it2=begin;
    it2++;
    while(it2!=end) {
      GradedMatrix& M1=*it1;
      GradedMatrix& M2=*it2;
      assert(M1.num_rows==M2.get_num_cols());
      M1.row_grades.clear();
      for(index j=0;j<M1.num_rows;j++) {
	M1.row_grades.push_back(M2.grades[j]);
      }
      //std::cout << "Copied " << M1.num_rows << " grades" << std::endl;
      it1++;
      it2++;
    }
    
  }



  template<typename InputIterator>
    void assign_grade_indices(InputIterator begin,
			      InputIterator end) {
    
    typedef typename InputIterator::value_type GradedMatrix;
    typedef typename GradedMatrix::Coordinate Coordinate;

    //std::cout << "Assgin grade indices with " << n1 << ", " << n2 << " columns" << std::endl;
    
    std::unordered_map<Coordinate,index> val_to_index_x, val_to_index_y;
    
    std::vector<Coordinate> x_vals,y_vals;

    for(auto it=begin;it!=end;it++) {
      GradedMatrix& matrix=*it;

      int n = matrix.get_num_cols();
      for(int i=0;i<n;i++) {
	x_vals.push_back(matrix.grades[i].first_val);
	y_vals.push_back(matrix.grades[i].second_val);
      }
      
    }
    
    std::sort(x_vals.begin(),x_vals.end());
    auto last_x = std::unique(x_vals.begin(),x_vals.end());
    x_vals.erase(last_x,x_vals.end());
    std::sort(y_vals.begin(),y_vals.end());
    auto last_y = std::unique(y_vals.begin(),y_vals.end());
    y_vals.erase(last_y,y_vals.end());
    
    if(verbose) std::cout << "Found " << x_vals.size() << " different x-values and " << y_vals.size() << " different y-values" << std::endl;
    for(auto it=begin;it!=end;it++) {
      GradedMatrix& matrix=*it;
      matrix.x_vals.clear();
      std::copy(x_vals.begin(),x_vals.end(),std::back_inserter(matrix.x_vals));
      matrix.y_vals.clear();
      std::copy(y_vals.begin(),y_vals.end(),std::back_inserter(matrix.y_vals));
      matrix.num_grades_x = x_vals.size();
      matrix.num_grades_y = y_vals.size();
      //std::cout << "num grade x: " <<  matrix.num_grades_x << std::endl;
      //std::cout << "num grade y: " <<  matrix.num_grades_y << std::endl;
    }
    
    for(index i=0;i<x_vals.size();i++) {
      val_to_index_x[x_vals[i]]=i;
    }
    for(index i=0;i<y_vals.size();i++) {
      val_to_index_y[y_vals[i]]=i;
    }
    long total_count=0;
    for(auto it=begin;it!=end;it++) {
      GradedMatrix& matrix = *it;
      int n = matrix.get_num_cols();
      //std::cout << "Handling " << n << " columns " << std::endl;
      total_count+=n;
      for(int i=0;i<n;i++) {
	matrix.grades[i].first_index=val_to_index_x[matrix.grades[i].first_val];
	assert(matrix.grades[i].first_index < matrix.num_grades_x);
	matrix.grades[i].second_index=val_to_index_y[matrix.grades[i].second_val];
	assert(matrix.grades[i].second_index < matrix.num_grades_y);
      }
    }
    //std::cout << "n1=" << n1 << std::endl;
    //std::cout << "n2=" << n2 << std::endl;
    if(verbose) std::cout << "N is " << total_count << std::endl;
    copy_grades_to_rows(begin,end);
    for(auto it=begin;it!=end;it++) {
      it->grade_indices_assigned=true;
    }
  }


  
  template<typename GradedMatrix>
    void assign_grade_indices_of_pair(GradedMatrix& M1, GradedMatrix& M2) {
    
    typedef typename GradedMatrix::Coordinate Coordinate;

    int n1=M1.get_num_cols();
    int n2=M2.get_num_cols();

    //std::cout << "Assgin grade indices with " << n1 << ", " << n2 << " columns" << std::endl;
    if(n1==0 && n2==0) {
      return;
    }
    
    std::unordered_map<Coordinate,index> val_to_index_x, val_to_index_y;
    
    std::vector<Coordinate> x_vals,y_vals;
    
    for(int i=0;i<n1;i++) {
      x_vals.push_back(M1.grades[i].first_val);
      y_vals.push_back(M1.grades[i].second_val);
    }
    for(int i=0;i<n2;i++) {
      x_vals.push_back(M2.grades[i].first_val);
      y_vals.push_back(M2.grades[i].second_val);
    }
    
    std::sort(x_vals.begin(),x_vals.end());
    auto last_x = std::unique(x_vals.begin(),x_vals.end());
    x_vals.erase(last_x,x_vals.end());
    std::sort(y_vals.begin(),y_vals.end());
    auto last_y = std::unique(y_vals.begin(),y_vals.end());
    y_vals.erase(last_y,y_vals.end());
    
    if(verbose) std::cout << "Found " << x_vals.size() << " different x-values and " << y_vals.size() << " different y-values" << std::endl;
    M1.x_vals=x_vals;
    M1.y_vals=y_vals;
    M1.num_grades_x = x_vals.size();
    M1.num_grades_y = y_vals.size();
    M2.x_vals=x_vals;
    M2.y_vals=y_vals;
    M2.num_grades_x = x_vals.size();
    M2.num_grades_y = y_vals.size();
    
    for(index i=0;i<x_vals.size();i++) {
      val_to_index_x[x_vals[i]]=i;
    }
    for(index i=0;i<y_vals.size();i++) {
      val_to_index_y[y_vals[i]]=i;
    }
    
    for(int i=0;i<n1;i++) {
      M1.grades[i].first_index=val_to_index_x[M1.grades[i].first_val];
      M1.grades[i].second_index=val_to_index_y[M1.grades[i].second_val];
      }
    for(int i=0;i<n2;i++) {
      M2.grades[i].first_index=val_to_index_x[M2.grades[i].first_val];
      M2.grades[i].second_index=val_to_index_y[M2.grades[i].second_val];
    }
    M1.row_grades.clear();
    for(int i=0;i<n2;i++) {
      M1.row_grades.push_back(M2.grades[i]);
    }
    M1.grade_indices_assigned=true;
    M2.grade_indices_assigned=true;
    
    //std::cout << "n1=" << n1 << std::endl;
    //std::cout << "n2=" << n2 << std::endl;
    if(verbose) std::cout << "N is " << n1+n2 << std::endl;
  }

  
  template<typename GradedMatrix>
    void assign_grade_indices(GradedMatrix& M1) {
    GradedMatrix M2;
    M2.set_dimensions(0,M1.num_rows);
    M2.grades=M1.row_grades;
    assign_grade_indices_of_pair(M1,M2);
  }
  
  template<typename GradedMatrix,typename OutStream>
    void print_in_scc_format(std::vector<GradedMatrix>& matrices, OutStream& outstr, bool extend_by_zero_matrix=false) {
    
    outstr << "scc2020\n2\n";
    for(index d=0;d<matrices.size();d++) {
      outstr << matrices[d].get_num_cols() << " ";
    }
    outstr << matrices[matrices.size()-1].num_rows << " ";
    if(extend_by_zero_matrix) {
	  outstr << "0";
    }
    outstr << "\n";

    for(index d=0;d<matrices.size()-1;d++) {
      GradedMatrix& M = matrices[d];
      M.print_in_rivet_format(outstr,false,false); // false for "print no header", false for "print no rows"
    }
    matrices.back().print_in_rivet_format(outstr,false,extend_by_zero_matrix);
      
  }
  
  // Re-sorts row and columns of the matrix to have co-lex order
  template<typename GradedMatrix>
    void to_colex_order(GradedMatrix& M,bool reorder_rows=true) {

    typedef typename GradedMatrix::Grade Grade;
    typedef Pre_column_struct<Grade> Pre_column;

    Sort_pre_column<Pre_column,typename GradedMatrix::Coordinate_traits::Compare> sort_pre_column;

    if(reorder_rows) {
      // Hack solution: Re-use the Pre_column_struct for this task, even though we have rows here
      
      std::vector<Pre_column> pre_rows;
      std::vector<index> empty_col;
      for(int i=0;i<M.row_grades.size();i++) {
	pre_rows.push_back(Pre_column(i,M.row_grades[i],empty_col));
      }
      std::sort(pre_rows.begin(),pre_rows.end(),sort_pre_column);
      std::vector<index> row_map;
      row_map.resize(M.num_rows);
      for(int i=0;i<pre_rows.size();i++) {
	row_map[pre_rows[i].idx]=i;
      }
      for(int i=0;i<M.row_grades.size();i++) {
	M.row_grades[i]=pre_rows[i].grade;
      }
      // Now iterate through the matrix and update all entries according to row_map
      for(int i=0;i<M.get_num_cols();i++) {
	std::vector<index> col;
	M.get_col(i,col);
	for(int j=0;j<col.size();j++) {
	  col[j]=row_map[col[j]];
	}
	std::sort(col.begin(),col.end());
	M.set_col(i,col);
      }
    }
    
    { // Now the columns
      std::vector<Pre_column> pre_cols;
      for(int i=0;i<M.grades.size();i++) {
	std::vector<index> col;
	M.get_col(i,col);
	pre_cols.push_back(Pre_column(i,M.grades[i],col));
      }
      std::sort(pre_cols.begin(),pre_cols.end(),sort_pre_column);
      for(int i=0;i<M.grades.size();i++) {
	M.grades[i]=pre_cols[i].grade;
	M.set_col(i,pre_cols[i].boundary);
      }
		  
    }
    assign_grade_indices(M);

  }



  
}//of namespace mpp_utils
