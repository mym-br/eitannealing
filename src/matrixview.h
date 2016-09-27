/*
 * matrixview.h
 *
 *  Created on: Set 27, 2016
 *      Author: thiago
 */
 
#ifndef _MATRIXVIEW_H_
#define _MATRIXVIEW_H_

#endif	// _MATRIXVIEW_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <QAbstractTableModel>

namespace internal {
    template <class MatrixType> class matrixViewModel : public QAbstractTableModel {
	    private:
		    const Eigen::EigenBase<MatrixType> &innermatrix;
	    public:
		    matrixViewModel(const Eigen::EigenBase<MatrixType> &m):innermatrix(m) {}
		    int rowCount( const QModelIndex & parent = QModelIndex() ) const {
			return this->innermatrix.rows();
		    }
		    int columnCount( const QModelIndex & parent = QModelIndex() ) const {
			return this->innermatrix.cols();
		    }
		    QVariant  data( const QModelIndex & index, int role = Qt::DisplayRole ) const {
			if (!index.isValid() || role != Qt::DisplayRole) return QVariant();
			return this->innermatrix.derived().coeff(index.row(), index.column());
		    }
		    QVariant  headerData ( int section, Qt::Orientation  orientation, int role = Qt::DisplayRole ) const {
			if (role == Qt::DisplayRole) return QString::number(section);
			return QVariant();
		    }
    };

    template <class MatrixType, unsigned int _Mode> class SparseSymetricMatrixViewModel : public QAbstractTableModel {
	    private:
		    typename MatrixType::Nested innerref;
	    public:
		    SparseSymetricMatrixViewModel(const Eigen::SparseSelfAdjointView< MatrixType, _Mode > &m):innerref(m.matrix()) {}
		    int rowCount( const QModelIndex & parent = QModelIndex() ) const {
			return this->innerref.rows();
		    }
		    int columnCount( const QModelIndex & parent = QModelIndex() ) const {
			return this->innerref.cols();
		    }
		    QVariant  data( const QModelIndex & index, int role = Qt::DisplayRole ) const {
			if (!index.isValid() || role != Qt::DisplayRole) return QVariant();
			if((_Mode & Eigen::UnitDiag) && (index.row()==index.column())) return 1.0f;
			if(((_Mode & Eigen::Lower) && (index.row()<index.column())) ||
			   ((_Mode & Eigen::Upper) && (index.row()>index.column())))
			  return this->innerref.coeff(index.column(), index.row());
			return this->innerref.coeff(index.row(), index.column());
		    }
		    QVariant  headerData ( int section, Qt::Orientation  orientation, int role = Qt::DisplayRole ) const {
			if (role == Qt::DisplayRole) return QString::number(section);
			return QVariant();
		    }
    };
}

template <class MatrixType, unsigned int mode = 0> QAbstractTableModel *makeMatrixTableModel(const Eigen::EigenBase<MatrixType> &m)
{
    return new internal::matrixViewModel<MatrixType>(m);
}


template <class MatrixType, unsigned int mode = 0> QAbstractTableModel *makeMatrixTableModel(const Eigen::EigenBase<Eigen::SparseSelfAdjointView<MatrixType, mode> > &m)
{
    return new internal::SparseSymetricMatrixViewModel<MatrixType, mode>(m.derived());
}

