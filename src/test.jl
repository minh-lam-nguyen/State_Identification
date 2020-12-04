include("data_structures/fiboheap.jl")

using .FibonacciHeaps

function test_fibo()
    h1 = MakeHeap()
    println("Made heap 1:"), print(h1)
    Insert(h1, "now")
    println("Inserted the word now into heap 1"), print(h1)
    h2 = MakeHeap()
    println("Made another heap 2.")
    t = Insert(h2, "time")
    println("Inserted the word time into heap 2:"), print(h2)
    h3 = Union(h1, h2)
    println("Made heap 3, union of heap 1 and heap 2:"), print(h3)
    println("The minimum of h3 is now \"$(Minimum(h3).value)\".")
    xkeys = [Insert(h3, x) for x in  ["all", "good", "men"]]
    println("Inserted 3 more into h3:"), print(h3)
    println("The minimum of h3 is now \"$(Minimum(h3).value)\".")
    println("The extracted min from heap 3 is: ", ExtractMin(h3).value)
    println("h3 is now:"), print(h3)
    println("Decrease key of heap 3 value $(xkeys[3].value) with the word come:")
    DecreaseKey(h3, xkeys[3], "come")
    print(h3)
    println("Delete node with value $(xkeys[3].value) from heap 3:")
    Delete(h3, xkeys[3])
    print(h3)
    println("The minimum of h3 is now: ", Minimum(h3).value)
end
