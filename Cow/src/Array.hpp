#ifndef CowArray_hpp
#define CowArray_hpp

#include <cstdlib>
#include <array>
#include <vector>
#include <functional>




namespace Cow
{
    class Array;
    class HeapAllocation;
    class RegionIterator;




    /**
    A class to manage large allocations on the heap using RAII standard.
    */
    class HeapAllocation
    {
    public:

        /**
        Create a null allocation.
        */
        HeapAllocation();

        /**
        Destructor.
        */
        ~HeapAllocation();

        /**
        Allocate a heap block of the given size. Bytes are *not* zero-
        initialized.
        */
        HeapAllocation (std::size_t numberOfBytes);

        /**
        Create a heap allocation from a std::string.
        */
        HeapAllocation (std::string content);

        /**
        Construct this memory block from a deep copy of another one.
        */
        HeapAllocation (const HeapAllocation& other);

        /**
        Move constructor.
        */
        HeapAllocation (HeapAllocation&& other);

        /**
        Assign this block the contents of another (deep copy).
        */
        HeapAllocation& operator= (const HeapAllocation& other);

        /**
        Move-assign this block the contents of another (steal data from other).
        */
        HeapAllocation& operator= (HeapAllocation&& other);

        /**
        Return the number of bytes in use.
        */
        std::size_t size() const;

        /**
        Return the contents of the buffer as a string.
        */
        std::string toString() const;

        /**
        Return a HeapAllocation whose binary data has opposite endian value as
        this one.
        */
        HeapAllocation swapBytes (std::size_t bytesPerEntry) const;

        void* begin() { return allocation; }

        const void* begin() const { return allocation; }

        template <class T> T& getElement (std::size_t index)
        {
            return static_cast<T*>(allocation)[index];
        }

        template <class T> const T& getElement (std::size_t index) const
        {
            return static_cast<T*>(allocation)[index];
        }

        template <class T> T* begin()
        {
            return static_cast<T*>(allocation);
        }

        template <class T> T* end()
        {
            return static_cast<T*>(allocation) + numberOfBytes / sizeof(T);
        }

        template <class T> const T* end() const
        {
            return static_cast<const T*>(allocation) + numberOfBytes / sizeof(T);
        }

    private:
        void* allocation;
        std::size_t numberOfBytes;
    };




    /**
    A type to represent the shape of an array or region.
    */
    using Shape = std::array<int, 5>;




    /**
    A type to represent a multi-dimensional index (i, j, k, m, n).
    */
    using Index = std::array<int, 5>;




    /**
    A helper class for manipulating a particular convention for 3D Array
    shapes. Axes 0, 1, and 2 are spatial, axis 3 is for scalar or vector
    components, and axis 4 (whose size is caled 'rank') is used to refer
    to indexed mesh locations (e.g. faces or edges) in a structured mesh.
    */
    class Shape3D
    {
    public:
        Shape3D();
        Shape3D (int n1, int n2, int n3);
        Shape3D (Shape S);
        Shape3D (const Array& A);
        operator Shape() const;
        const int &operator[] (int index) const;
        int &operator[] (int index);
        Shape3D operator*(int x) const;
        Shape3D operator/(int x) const;
        /** Shape with spatial axes reduced. */
        Shape3D reduced (int delta=1) const;
        /** Shape with spatial axes reduced by first 3 elements of delta. */
        Shape3D reduced (Shape delta) const;
        /** Shape with the given axis reduced by delta. */
        Shape3D reduced (int axis, int delta) const;
        /** Shape with the given axis increased by delta. */
        Shape3D increased (Shape delta) const;
        /** Shape with spatial axes increased by first 3 elements of delta. */
        Shape3D increased (int delta=1) const;
        /** Shape with the given axis increased by delta. */
        Shape3D increased (int axis, int delta) const;
        /** Shape with axis 3 replaced. */
        Shape3D withComponents (int numComponents) const;
        /** Shape with axis 4 replaced. */
        Shape3D withRank (int rank) const;
        /** True if other is <= in size to this on each axis. */
        bool contains (Shape3D other) const;

        /**
        A utility function which deploys a function of i, j, k over the given
        shape. This is essentially short-hand for writing a triple for-loop.
        */
        void deploy (std::function<void (int i, int j, int k)> function) const;

    private:
        Shape S;
    };




    /**
    A class that represents a relative range along a single array axis.
    */
    class Range
    {
    public:
        /**
        Construct a relative or absolute range [i0:i1:di]. The default [0:0:1]
        is a relative range that covers the whole extent.
        */
        Range (int lower, int upper, int stride=1);

        /**
        Construct a relative range from the character ":" as a shortcut for
        the whole extent.
        */
        Range (const char*);

        /** Return true if the upper bound is relative to end. */
        bool isRelative() const;

        /**
        Get the extent of an absolute range, this is simply upper - lower.
        */
        int extent() const { return upper - lower; }

        /**
        Get the size of an absolute range, after strides are accounted for.
        */
        int size() const;

        /**
        Get the absolute size of a relative range for the given axis size,
        after strides are accounted for.
        */
        int size (int absoluteSize) const;

        /** Return an aboslute version of this range. */
        Range absolute (int absoluteSize) const;

        const int lower;
        const int upper;
        const int stride;
    };


    /**
    A class that represents a relative or absolute region within an array. If
    any component of upper is either zero or negative, then that value
    indicates a distance from the end of the array, and the region is called
    'relative'. An 'absolute' region is generated by providing a shape object
    to the absolute() method. By default each axis covers the range [0:0:1],
    which denotes the entire extent of the axis. [0:0:0] is reserved for the
    empty region.
    */
    class Region
    {
    public:

        /**
        Return an object that denotes an empty region.
        */
        static Region empty();

        /**
        Return an absolute region for the given shape.
        */
        static Region whole (Shape shape);

        /**
        Construct a default region, which is relative and refers to the entire
        extent of its target. Alias for "whole".
        */
        Region (Shape shape);

        /**
        Construct a default region, which is relative and refers to the entire
        extent of its target.
        */
        Region();

        /**
        Return a new region, with the upper, lower, and stride on the given
        axis replaced.
        */
        Region withRange (int axis, int lower, int upper, int stride=1) const;

        /**
        Return a new region, with the lower bound on the given axis replaced.
        */
        Region withLower (int axis, int newLower) const;

        /**
        Return a new region, with the upper bound on the given axis replaced.
        */
        Region withUpper (int axis, int newUpper) const;

        /**
        Return a strided version of this region, along the given axis.
        */
        Region withStride (int axis, int newStride) const;

        /**
        Return true if the upper bound is relative to end.
        */
        bool isRelative() const;

        /**
        Check if this region is empty.
        */
        bool isEmpty() const;

        /**
        Check if two regions are identical.
        */
        bool operator== (const Region& other) const;

        /**
        Return the total number of elements in the region, after strides
        accounted for. The region is assumed to be absolute.
        */
        int size() const;

        /**
        Return the number of elements along each axis, after strides are
        accounted for. The region is assumed to be absolute.
        */
        Shape shape() const;

        /**
        Return a 3D version of shape().
        */
        Shape3D shape3D() const;

        /**
        Return shape(), but with trailing axes of length 1 removed.
        */
        std::vector<int> getShapeVector() const;

        /**
        Return an absolute version of this region by providing a definite shape.
        */
        Region absolute (Shape shape) const;

        /**
        Return an absolute version of this region by providing a shape vector.
        Trailing axes are assumed to have size equal to 1.
        */
        Region absolute (std::vector<int> shapeVector) const;

        /**
        If necessary, modify this region in-place to be absolute.
        */
        void ensureAbsolute (Shape shape);

        /**
        Return the absolute range of indices (with stride information) covered
        for the given axis.
        */
        Range range (int axis) const;

        Index lower;
        Index upper;
        Index stride;
    };




    /**
    A multidimensional array class, hard-coded to accommodate up to 5 axes.
    */
    class Array
    {
    public:
        class Reference;
        class Iterator;

        /**
        Various array constructors. Arrays are zero-initialized on construciton.
        */
        Array();
        Array (Shape shape);
        Array (Reference reference);
        Array (int n1);
        Array (int n1, int n2);
        Array (int n1, int n2, int n3);
        Array (int n1, int n2, int n3, int n4);
        Array (int n1, int n2, int n3, int n4, int n5);

        /**
        Copy constructor.
        */
        Array (const Array& other);

        /**
        Move constructor.
        */
        Array (Array&& other);

        /**
        Assignment operator.
        */
        Array& operator= (const Array& other);

        /**
        Move-assignment operator.
        */
        Array& operator= (Array&& other);

        /** Get a reference to the underlying buffer. */
        HeapAllocation& getAllocation() { return memory; }

        /** Get a const reference to the underlying buffer. */
        const HeapAllocation& getAllocation() const { return memory; }

        /**
        Return the total number of doubles in this array.
        */
        int size() const;

        /**
        Return the number of elements along one of the axes.
        */
        int size (int axis) const;

        /**
        Return the Array's shape as a 5-component array.
        */
        Shape shape() const;

        /**
        Return a 3D version of shape().
        */
        Shape3D shape3D() const;

        /**
        Return the strides in memory along each axis.
        */
        Shape strides() const;

        /**
        Return the shape of this array as a vector, with trailing axes that
        have length 1 removed.
        */
        std::vector<int> getShapeVector() const;

        /**
        Return a new array that is the transpose of this one, A.transpose()
        (i, j, k, m, n) == A (n, m, k, j, i). The returned array has same
        memory layout type as this.
        */
        Array transpose() const;

        /**
        Return an array with the pair of given axes transposed.
        */
        Array transpose (int axis1, int axis2) const;

        /**
        Extract a deep copy of the given relative or absolute region of this
        array. This is equivalent to auto B = Array (A[region]);
        */
        Array extract (Region R) const;

        /**
        Insert all of the source array into the given region of this array.
        The region argument may be omitted, in which case this array's data is
        replaced completely. This operation may improve efficiency relative to
        the assignment operator, since this array's HeapAllocation is reused.
        */
        void insert (const Array& source, Region R=Region());

        /**
        Copy data from the source region of A into the target region of this
        array.
        */
        void copyFrom (const Array& A, Region targetRegion, Region sourceRegion=Region());

        /**
        Change the Array's shape, without modifying its data layout. The new
        size must equal the old size.
        */
        void reshape (int n1, int n2=1, int n3=1, int n4=1, int n5=1);
        
        /**
        Return a trivial iterator to the beginning of the array.
        */
        double* begin() { return memory.begin<double>(); }

        /**
        Return a trivial iterator to the end of the array.
        */
        double* end() { return memory.end<double>(); }

        /**
        Return a trivial iterator to the end of the array.
        */
        const double* end() const { return memory.end<double>(); }

        /** Retrieve a value by linear index. */
        double& operator[] (int index);

        /** Retrieve a const value by linear index. */
        const double& operator[] (int index) const;

        /**
        Return a reference to a particular region in this array. It is the
        caller's responsibility to ensure the referenced array remains alive
        at least as long as the reference does.
        */
        Reference operator[] (Region R);

        operator Reference() { return this->operator[] (Region()); }

        double& operator() (int i);
        double& operator() (int i, int j);
        double& operator() (int i, int j, int k);
        double& operator() (int i, int j, int k, int m);
        double& operator() (int i, int j, int k, int m, int n);

        const double& operator() (int i) const;
        const double& operator() (int i, int j) const;
        const double& operator() (int i, int j, int k) const;
        const double& operator() (int i, int j, int k, int m) const;
        const double& operator() (int i, int j, int k, int m, int n) const;

        /**
        Return an array which results from applying the given callback
        function element-wise.
        */
        Array map (std::function<double (double)> function) const;

        static Shape shapeFromVector (std::vector<int> shapeVector);

        static std::vector<int> vectorFromShape (Shape shape);

        static bool isBoundsCheckDisabled();
        
        /**
        A utility function which deploys a function over the first three axes
        of the given shape. This is essentially short-hand for writing a
        triple for-loop.

        \deprecated
        move to Shape3D class
        */
        static void deploy (Shape shape, std::function<void (int i, int j, int k)> function);

        class Reference
        {
        public:
            /**
            Constructor. The region is assumed to be absolute.
            */
            Reference (Array& A, Region R);

            /**
            Copy values from a source array to the referenced region of the
            target array.
            */
            const Array& operator= (const Array& source);

            /**
            Copy from a reference into another array.
            */
            const Array::Reference& operator= (const Array::Reference& source);

            /**
            Return the referenced array.
            */
            Array& getArray();
            const Array& getArray() const;

            /**
            Return the region of the referenced array.
            */
            const Region& getRegion() const;

            /**
            Convenience for shape()[axis].
            */
            int size (int axis) const;

            /**
            Return the referenced regions's shape, short for ref.getRegion().shape().
            */
            Shape shape() const;

            /**
            Convenience function for getRegion().getShapeVector().
            */
            std::vector<int> getShapeVector() const;

            /**
            Return an iterator to the beginning of the array. When
            incremented, the result is equivalent to a multidimensional loop
            over the referenced rgion.
            */
            Iterator begin();

            /**
            End iterator.
            */
            Iterator end();
        private:
            friend class Array;
            Array& A;
            Region R;
        };

        class Iterator
        {
        public:
            Iterator (Array& A, Region R, bool isEnd=false);

            /**
            Treat the iterator as a double pointer to the value at the current index.
            */
            operator double*() const;

            /** Increment operator. */
            double* operator++ ();

            /**
            Return the value at some distance in memory away form the current
            location.
            */
            double& operator[] (int offset);

            /** Comparison operator. */
            bool operator== (const Iterator& other) const;

            /** Print the current index of the iterator. */
            void print (std::ostream& stream) const;

            /**
            Get the current index, relative to the beginnning of the
            referenced array.
            */
            Index index() const;

            /**
            Return the current index of the iterator relative to the start of
            the referenced region, rather than the underlying array.
            */
            Index relativeIndex() const;

        private:
            double* getAddress() const;
            Array& A;
            Region R;
            Index currentIndex;
            double* currentAddress;
        };

    private:
        /** @internal */
        static void copyRegion (Array& dst, const Array& src, Region source, Region target);
        int n1, n2, n3, n4, n5;
        Shape S;
        HeapAllocation memory;
    };
};


std::ostream& operator<< (std::ostream &stream, const Cow::HeapAllocation &memory);

#endif