class Test1:
    def __init__(self):
        self.attribute1 = "Hello world!"

class Test2(Test1):
    def __init__(self, Test1):
        self.attribute1 = Test1.attribute1
        self.attribute2 = "World, hello!"


obj1 = Test1()
obj2 = Test2(obj1)
print(obj1.attribute1)
