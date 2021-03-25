    def decode_json(self, file_, decoder=None,
                deserialize_type=None):
        decoder = str(decoder).lower()
        deserialize_type = str(deserialize_type).lower()
        decode_packages = ['json', 'sjson']
        deserializers = ['load', 'loads']
        assert decoder in decode_packages, \
            f"{decoder} not in \{decode_packages}"
        assert deserialize_type in deserializers, \
        f"{deserialize_type} not in {deserializers}"
        if decoder == 'sjson':
            if deserialize_type == 'load':
                data = sjson.load(file_)
            elif deserialize_type == 'loads':
                data = sjson.loads(file_)
        elif decoder == 'json':
            if deserialize_type == 'load':
                data = json.load(file_)
            elif deserialize_type == 'loads':
                data = json.loads(file_)
        if data is None:
            print("No data read")
        return data