#ifndef OBJECT_POOL
#define OBJECT_POOL

template <typename T>
class Object_pool {
private:
    std::deque<T*> pool;
public:
    template <typename... Args>
    T* acquire(Args&&... args) {
        if (pool.empty()) {
            return new T(std::forward<Args>(args)...);
        }
        T* obj = pool.back();
        pool.pop_back();
        obj->reset(std::forward<Args>(args)...); // Initialize the object
        return obj;
    }

    void release(T* obj) {
        pool.push_back(obj);
    }

    ~Object_pool() {
        for (T* obj : pool) {
            delete obj;
        }
    }
};

#endif // !OBJECT_POOL
